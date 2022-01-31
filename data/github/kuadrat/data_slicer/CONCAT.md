![](https://raw.githubusercontent.com/kuadrat/data_slicer/master/img/02_Data_Slicer_Logo.png)

# Data slicer 
The `data-slicer` package offers fast tools for inspection, visualization, 
slicing and analysis of 3(+) dimensional datasets at a general level.
It also provides a framework and building blocks for users to easily create 
more specialized and customized tools for their individual use cases.

![](https://raw.githubusercontent.com/kuadrat/data_slicer/master/img/pit_demo.gif)

`data-slicer` was originally developed to deal with the high data throughput of 
modern measurement instruments, where quick visualizations and preliminary 
analyses are necessary to guide the direction of a measurement session.
However, the package is designed to be agnostic of the concrete use-case and 
all scientific, engineering, medical, artistic or other data driven 
disciplines where inspection and slicing of (hyper)cubes is required could 
potentially benefit from `data-slicer`.

## Documentation

This README just gives a minimal overview.
For more information, guides, examples and more, visit the documentation which is hosted 
by the friendly people over at ReadTheDocs.org:
https://data-slicer.readthedocs.io/en/latest/

## Installation

`data-slicer` should run on all platforms that support python and has been 
shown to run on Windows, macOS ans Linux.

The package can be installed from [PyPI](https://pypi.org/project/data-slicer/) 
using `pip install data_slicer`.
It is recommended to do this from within some sort of virtual environment.
Visit the documentation for more detailed instructions:
https://data-slicer.readthedocs.io/en/latest/installation.html

### Dependencies

This software is built upon on a number of other open-source frameworks.
The complete list of packages can be found in the file `requirements.txt`.
Most notably, this includes [matplotlib](https://matplotlib.org/), 
[numpy](https://numpy.org/) and 
[pyqtgraph](https://pyqtgraph.readthedocs.io/en/latest/).

## Citing

If you use `data-slicer` in your work, please credit it by citing the following publication:

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02969/status.svg)](https://doi.org/10.21105/joss.02969)

Kramer et al., (2021). Visualization of Multi-Dimensional Data -- The data-slicer Package. Journal of Open Source Software, 6(60), 2969, https://doi.org/10.21105/joss.02969

## Contributing

You are welcome to help making this software package more useful!
You can do this by giving feedback, reporting bugs, issuing feature requests 
or fixing bugs and adding new features yourself. 
Furthermore, you can create and share your own plugins (refer to the 
[documentation](https://data-slicer.readthedocs.io/en/latest/)).

### Feedback, bugs and feature requests

The most straightforward and organized way to help improve this software is 
by opening an **issue** on [the github 
repository](https://github.com/kuadrat/data-slicer/issues).
To do this, navigate to the **Issues** tab and click **New issue**.
Please try to describe the bug you encountered or the new feature you would 
like to see as detailed as possible.
If you have several bugs/ideas please open a separate issue for each of them 
in order to keep the discussions focused.

If you have anything to tell me that does not seem to warrant opening an 
issue or you simply prefer to contact me directly you can do this via e-mail:
kevin.pasqual.kramer@protonmail.ch

### Contributing code

If you have fixed a bug or created a new feature in the source code yourself, 
it can be merged into this project.
Code contributions will be acknowledged in this README or, if the number of 
contributors grows too large, in a separate file.
If you are familiar with the workflow on github, please go ahead and create a 
pull request.
If you are unsure you can always contact me via e-mail (see above).

### Plugins

If you have created a PIT plugin, feel free to add it to the list below, 
either via pull request or through an e-mail (see above).
Also check the 
[documentation](https://data-slicer.readthedocs.io/en/latest/contributing.html) 
for a guide on how to create a plugin.

| Plugin Name (link) | Description | 
| ------------------ | ----------- |
| [ds-example-plugin](https://github.com/kuadrat/ds-example-plugin) | exists as a minimal example that can be used for guidance when creating your own plugins. A step-by-step tutorial on how it was made can be found in [the documentation](https://data-slicer.readthedocs.io/en/latest/contributing.html) |
| [ds-arpes-plugin](https://github.com/kuadrat/ds_arpes_plugin) | tools for angle-resolved photoelectron spectroscopy (ARPES) | 

---
title: Visualization of Multi-Dimensional Data -- The data-slicer Package
bibliography: paper.bib
tags:
  - Python
  - visualization
authors:
  - name: Kevin Kramer
    orcid: 0000-0001-5523-6924
    affiliation: 1
  - name: Johan Chang
    orcid: 0000-0002-4655-1516
    affiliation: 1
affiliations:
  - name: Physik Institut, Universität Zürich, Winterthurerstrasse 190, CH-8057 Zürich, Switzerland
    index: 1
date: 25 February 2021
---

# Statement of Need
\label{sec:intro}

From prehistoric cave-wall paintings to the invention of print and most 
recently electronic hard-disks, human data storage capacity has evolved 
tremendously.
Information/data is of great value and hence associated with innovation and 
technological progress.
This is especially true in analytical disciplines i.e. all sciences ranging 
from physics to psychology and medicine.
In observational sciences, most measurement techniques undergo steady 
improvements in acquisition time and resolution.
As a result the sheer data throughput is continually increasing.
Examples of techniques where the typical data output has moved from 1D to 3D 
in the past few decades are shown in \autoref{fig1}.

More data is always welcome.
However, in many disciplines human digestion of these large amounts of data 
has now become the bottleneck.
In many fields, for example those working at large scale synchrotron 
facilities where the duration of the experiment is limited, scientists 
require a means of quick data inspection and carrying out a fast preliminary 
analysis in order to take decisions on the course of the experiment.
Many of the existing powerful and versatile visualization tools 
[@fedorov123d;@ahrens05paraview;@00mayavi;@00visit] are not well suited for 
this purpose and cannot easily support the specific and often changing needs 
of a given discipline.
The result is that each community ends up developing their own solutions to 
the problem of quick data visualization and inspection, e.g.
[@stansbury20pyarpes;@lass20mjolnir].
However, since these implementations are usually intertwined and 
entangled with the community-specific parts, such solutions are typically 
not transferrable across different disciplines or experimental methodologies.

We have developed PIT and the data-slicer package to account 
for these needs: offering tools for fast live visualization of data at a 
general scope that can easily be adjusted and fine tuned for more 
specific problems.

![Evolution of data acquisition in the field of spectroscopy. 
(a,b) Angle resolved photoemission electron spectroscopy (ARPES) 
[@shai13quasiparticle; @wells92evidence], 
(c,d) tunnelling spectroscopy (STS) [@zhang19machine; @giaever62tunneling], 
and (e,f) inelastic neutron scattering (INS) 
[@Wan_2020; @Bastien; @woods60lattice] spectroscopy 
techniques all started with single spectrum collection (top row).
Modern spectroscopic and scattering techniques, however, involve 
multidimensional data acquisition (bottom row).
  \label{fig1}
](fig1.pdf)

# Summary

data-slicer is a python package that contains several functions and classes 
which provide modular Qt [@company00qt; @00riverbank] widgets, tools and 
utilities for the visualization of three-dimensional (3D) datasets.
These building blocks can be combined freely to create new applications.
Some of these building blocks are used within the package to form a 
graphical user interface (GUI) for 3D data visualization and manipulation: 
the Python Image Tool (PIT).
The relation between different elements of the package and external software 
is schematically depicted in \autoref{fig2}.

## PIT
PIT consists of a number of dynamic plot figures which allow browsing through 
3D data by quickly selecting slices of variable thickness from the data cube 
and further cutting them up arbitrarily.
Two core features of PIT should be explicitly mentioned.
The first is the ability to create slices of the 3D data cube along arbitrary 
angles quickly.
This is facilitated on the GUI side through a simple draggable line to select 
the slice direction.
The superior speed of this operation is enabled by the use of optimized 
functions.
The second feature worth mentioning is the inclusion of an ipython console 
which is aware of the loaded data as well as of all GUI elements.
The console immediately enables users to run any analysis routine suitable to 
their respective needs.
This includes running python commands in a workflow familiar to pylab or 
Jupyter [@00jupyter] notebook users but also loading or directly running 
scripts into or from the console, using ipython's line magic functions 
\texttt{\%load} and \texttt{\%run} respectively.
Effectively, this design is central in empowering users to accomplish any task
imaginable --- as long as it is possible to achieve with python.

## Plugins
It is clear that it can get complicated and tedious to run certain types of 
data processing or analysis from the ipython console, as described in the 
previous paragraph.
For such cases, PIT provides an additional level of customizability and 
control through its plugin system.
Plugins are regular python packages that can be loaded from within PIT and 
enhance it with new functionality.
A plugin can interact with all elements in PIT via the same interfaces as can 
be done through the built-in ipython console.
Creating a plugin therefore requires little programming skills and no further 
knowledge of the inner workings of PIT.
In this manner, different communities of users can create and share their 
field-specific plugins which allow them to customize PIT to their needs.

As an example, we mention the ds-arpes-plugin [@arpesPlugin], which provides 
basic functionalities for loading of ARPES datasets and handles for typical 
analysis functions, customized and taylored to be used from within PIT.

## Modularity and widgets

PIT is constructed in a modular fashion, constituting of different widgets 
that have been combined together to make a useful, ready to use tool.
However, different applications may require slightly different 
functionalities, and the setup in PIT may not be optimal.
The data-slicer package makes all the used widgets in PIT and some additional ones 
independently available to the user.
These widgets can be arbitrarily combined to create customized applications 
in a relatively simple manner.

In summary, the data-slicer package solves the problem of offering the right 
scope -- neither too specialized that it can only be used by a narrow community 
nor too bloated such that it becomes hard to do specific operations -- by 
offering a variety of methods for users of varying backgrounds to get exactly 
the tools they need.
On the first and most general level, PIT offers a ready-to-use GUI for quick 
3D data visualization without any need of programmatic user input.
Users can satisfy their more specific needs either through use of the console 
or by implementing a plugin, which can both be accomplished with little 
programming knowledge.
On the last, most specific level users can use and arrange the building 
blocks contained in the package to create completely new applications or 
embed PIT or other parts of the data-slicer package into an existing application.

![Schematic structural overview of the data-slicer package (represented by 
the large grey block). 
The python image tool (PIT) is made up of different modular building blocks. 
These blocks can be used to create or enhance new applications.
PIT itself can be used directly, optionally augmented through plugins.
Alternatively, it can be embedded in external applications.
\label{fig2}
](fig2.pdf)

# Acknowledgements

We are thankful for fruitful discussions with and inputs from Titus Neupert, 
Claude Monney, Daniel Mazzone, Nicholas Plumb, Wojciech Pude\l{}ko, Niels 
Bech Christensen as well as for the testing efforts of Qisi Wang, Julia 
Küspert and Karin von Arx.
We acknowledge the contribution of other open-source frameworks to this 
publication, namely matplotlib [@hunter07matplotlib], numpy 
[@harris20array] and pyqtgraph [@00pyqtgraph].
Kevin Kramer and Johan Chang acknowledge support by the Swiss National 
Science Foundation.


# References
.. _sec-contributing:

Contributing
============

Thank you for helping to make this software package more useful!
There are several ways in which you can contribute: giving feedback, 
reporting bugs, issuing feature requests or fixing bugs and adding new 
features yourself. 
Furthermore, you can create and share your own :ref:`plugins <sec-plugin>`.

Feedback, bugs and feature requests
-----------------------------------

The most straightforward and organized way to help improve this software is 
by opening an **issue** on `the github 
repository <https://github.com/kuadrat/data-slicer/issues>`_.
To do this, navigate to the **Issues** tab and click **New issue**.
Please try to describe the bug you encountered or the new feature you would 
like to see as detailed as possible.
If you have several bugs/ideas please open a separate issue for each of them 
in order to keep the discussions focused.

If you have anything to tell me that does not seem to warrant opening an 
issue or you simply prefer to contact me directly you can do this via e-mail:
kevin.pasqual.kramer@protonmail.ch

Contributing code
-----------------

If you have fixed a bug or created a new feature in the source code yourself, 
it can be merged into this project.
Code contributions will be acknowledged in the README or, if the number of 
contributors grows too large, in a separate file.
If you are familiar with the workflow on github, please go ahead and create a 
pull request.
If you are unsure you can always contact me via e-mail (see above).

Plugins
-------

If you have created a PIT plugin, feel free to add it to the list on the 
:ref:`plugins page <sec-plugin>`, either via pull request or through an 
e-mail (see above).

.. _sec-plugin-how-to:

Creating a plugin
.................

Creating a plugin is not difficult, but there are a few details that one 
needs to get right in order for things to work.
The package `ds-example-plugin <https://github.com/kuadrat/ds-example-plugin>`_
was made to serve as an example to showcase how a plugin can be created.
This section gives a step-by-step guide on how ``ds-example-plugin`` was 
built and should allow you to figure out how to create your own plugins by 
analogy.
It is helpful to be familiar with the concept of object oriented 
programming (working with classes and inheritance), but not mandatory.

Before we start, let's have a look at what the example plugin does.
You can either install it using ``pip install ds-example-plugin`` or download 
the source code from github and save it (or a link to it) in the ``plugins`` 
subdirectory of your :ref:`config directory <sec-config>`.
Once downloaded and installed, start PIT and load the plugin::

   [1] example = mw.load_plugin('ds_example_plugin')
   Importing plugin ds_example_plugin (Example plugin).
   
Running ``example.help()`` then gives us a list of available functions. 
There are just two (besides the ``help`` function).
Let's try them both::

   [2] example.example_function()
   This message is sent to you by the Example plugin.

   [3] example.blur()

The first one doesn't do anything but print some text to the console.
The second one blurs the data.
It can also be called with a numeric argument to specify the *amount* of 
blurring.
Since the blurring is applied directly to the data, it can only be undone by 
using ``pit.reset_data()``.

That's really all this plugin does.
Now let's see how it is implemented.
In the following, we refer to the source code files that are found `on the 
github <https://github.com/kuadrat/ds-example-plugin>`_.

Step 1: Writing the Plugin class
********************************

The file ``ds_example_plugin.py`` contains the core of our plugin.
Essentially, we define a ``class`` that inherits from 
:class:`data_slicer.plugin.Plugin`::

   from data_slicer import plugin

   class Example_Plugin(plugin.Plugin) :
      [...]

The ``__init__()`` method of this class **must** call the superclass' 
``__init__()`` through the line::

   def __init__(self, *args, **kwargs) :
      [...]
      super().__init__(*args, **kwargs)
      [...]

It is recommended that you give your plugin a ``name`` and ``shortname`` by 
setting the respective attributes in the ``__init__()`` function::

   def __init__(self, *args, **kwargs) :
      [...]
      self.name = 'Example plugin'
      self.shortname = 'example'

Of course you can add any kind of code you need to in the ``__init__()`` 
function.
But to summarize, the bare requirements that absolutely need to be fulfilled 
for a plugin to work are:

   1. Inheriting from :class:`data_slicer.plugin.Plugin`.

   2. Calling the superclass' constructor with ``super().__init__(...)``

Step 2: Adding functionality to the plugin class
************************************************

We can now add functions to this class.
As a first example, let's have a look at ``example_function``::

   def example_function(self) :
      """ [...] """
      print([...])

The function is defined like any python function, but as a class function it 
**must** contain the argument ``self``.

The *doc string* (everything that appears between the ``"""`` directly after 
the function name) will be visible to the user when they use the built-in 
``help()`` function, so it's a good idea to explain as clearly as possible 
what the function does and how to use it properly.

As we have seen, this function is later available to the user through 
``example.example_function()`` and it will also be listed by ``example.help()``.

However, this ``example_function()`` doesn't actually do anything useful and 
related to PIT, so let's have a look at ``blur()`` instead, where we see how 
PIT's data and other elements can be accessed::

   # Add to the imports on top
   from scipy.ndimage import gaussian_filter

   [...]

   # Inside the class definition
      def blur(self, sigma=1) :
         """ [...] """
         data = self.data_handler.get_data()
         self.data_handler.set_data(gaussian_filter(data, sigma))

We see again how the argument ``self`` needs to appear first.
After it, arbitrary positional and keyword argument can follow, as usual with 
python functions.

In the first line of the function definition we see that we can get access to 
the data through ``self.data_handler.get_data()``.
``data_handler`` is available to us thanks to inheriting from 
:class:`data_slicer.plugin.Plugin` and is the same as ``pit`` which is 
available from the ipython console, i.e. it is a 
:class:`data_slicer.pit.PITDataHandler` object.
Similarly, within our class definition we also have access to a 
:class:`data_slicer.pit.MainWindow` object through ``self.main_window``.
It is through these two objects that we can access all the data and visual 
elements of PIT.
Check out their respective documentations to see what functions they provide, 
but the ones that are likely most useful are 
:meth:`self.data_handler.get_data() <data_slicer.pit.PITDataHandler.get_data>`,
:meth:`self.data_handler.set_data() <data_slicer.pit.PITDataHandler.set_data>`
and :meth:`self.data_handler.overlay_model() 
<data_slicer.pit.PITDataHandler.overlay_model>`.

So to summarize step 2, we can add arbitrary functionality to our plugin by 
defining functions.
By means of ``self.data_handler`` and ``self.main_window`` we get access to 
elements of PIT, meaning we can read and change them.
Every function we define will be directly available to the user and gets 
automatically listed by the plugin's ``help()`` function.

Step 3: Creating the module structure
*************************************

PIT loads plugins by using python's ``import`` capabilities to import the 
plugin class we created in steps 1 and 2.
For this to work, python needs to recognize our code as a module.
To be recognized as a module, we simply need to have a file called 
``__init__.py`` in the same directory as our code.
For python it is enough for this file to exist, even if it is empty.
However, if you look at the ``__init__.py`` in the code of 
``ds-example-plugin`` you'll find that it does contain one line of code::

   from ds_example_plugin.ds_example_plugin import Example_Plugin as main

This line is necessary for PIT's plugin mechanism to work.
It simply imports the class we defined in steps 1 and 2 under the alias 
``main``.

Let's have a closer look at the structure of this line, which can be stripped 
down to consist of four parts::

   from <1>.<2> import <3> as <4>

When loading the plugin, PIT will look for a class called ``main`` and 
instantiate it - this means that when you write your plugin, the last word in 
this line (``<4>``) always has to be ``main``.

Part ``<3>`` must match the name of the class from earlier and part ``<2>`` 
is the exact name of the file (minus the ``.py`` suffix) that contains said 
class definition.
Finally, part ``<1>`` is the exact name of the directory that contains both 
the ``__init__.py`` file and the file ``<2>``.
``<1>`` is also the name that will have to be used to load the plugin from 
within PIT with ``mw.load_plugin(<1>)``.
``<1>`` and ``<2>`` don't necessarily have to be the same, but it is quite 
customary to use the same name here.

In summary, in this step we created the file ``__init__.py`` that contains 
just a very specific line of code.

Final step: make the plugin available
*************************************

In order for PIT to find your plugin, it has to

   - either be detectable by python

   - or be placed in your :ref:`config directory <sec-config>`.

To be detectable by python, you have to place the directory containing the 
plugin file and the ``__init__.py`` file (directory with name ``<1>`` from 
above) somewhere in your PYTHONPATH.
If you're just using a plugin for personal use, it is probably easiest to 
simply place it in your config directory.
Don't forget that you have the option of letting PIT :ref:`automatically load 
plugins <sec-config-plugins>`, using the ``autoload.txt`` file.

If you want to share your plugin with other people, it is enough to make the 
source code available by any means.
Other users can then just download the code and place it in their config 
directories to make use of the plugin.
However, if you know how to properly package a python module, you could do 
that and make it available on PyPI, such that users can simply install the 
plugin by getting it with ``pip``.
This would have the advantage that it will make it easier for users to get 
updated versions of the plugin.
A tutorial on packaging would, however, go beyond the scope of this tutorial.

Don't forget to contact me if you have created and shared a plugin, such that 
I can add it to the :ref:`list of available plugins <sec-plugin-list>`!

.. _sec-data_format:

Data formats and data loading
=============================

Accessing the displayed data
----------------------------

The data currently represented in ``PIT`` can always be accessed and changed 
from the ``ipython`` console through ``pit.get_data()`` and ``pit.set_data()``.

We could, for example, replace the first slice of our data by random numbers::

   # Get hold of the data
   [1] data = pit.get_data()
   # Get the shape
   [2] x, y, z = data.shape
   # Randomize the first slice
   [3] data[:,:,0] = np.random.rand(x, y)
   # Set the data again
   [4] pit.set_data(data)

The :meth:`~data_slicer.pit.PITDataHandler.set_data` function takes an 
additional argument **axes**, which should be a list/tuple/array of 3 arrays 
representing the x, y and z axis coordinates.
The following example would rescale the x-axis by an order of magnitude and 
shift the y-axis by 50 units::

   # Get the current axes
   [1] axes = pit.axes
   # Set modified axes coordinates with the set_data method
   [2] pit.set_data(axes=[axes[0]*10, axes[1]+50, axes[2]])

As another example, this is how we could apply some Gaussian blurring to the 
data (requires ``scipy`` to be installed)::

   # Import the Gaussian filter from scipy
   [1] from scipy.ndimage import gaussian_filter
   # Get the current data
   [2] data = pit.get_data()
   # Set blurred data
   [3] pit.set_data(gaussian_filter(data, 1))
   # Play with different levels of blurring
   [4] pit.set_data(gaussian_filter(data, 10))

If you have certain operations that you routinely carry out on your data, it 
is recommended to automate this process by writing a :ref:`plugin 
<sec-plugin>`.

Loading data from a file
------------------------

The existing dataloaders found in the :mod:`dataloading module 
<data_slicer.dataloading>` can handle two types of input data formats:

1. A binary file that has been created using python's ``pickle`` module.
   The pickled object can be either a ``numpy`` array, or a dictionary or 
   :class:`argparse.Namespace` containing the keys *data* and *axes*, as 
   described in :class:`~data_slicer.dataloading.Dataloader_Pickle`.

2. A plain text file containing values for x, y, z and the actual data in 
   four columns as described in the documentation of 
   :class:`~data_slicer.dataloading.Dataloader_3dtxt`.
   Notice that there is an utility function 
   (:func:`~data_slicer.dataloading.three_d_to_txt`) that can help you 
   create such a text file in the correct format from existing data.

In the following we give an example tutorial for both cases.


Binary pickle files
'''''''''''''''''''''''''''''''''

Step 1: Load the data into python
.................................

We will assume our starting point to be a 3D data file of any format.
It could be an `Igor file <https://www.wavemetrics.com/igor-8-highlights>`_, 
an `HDF5 file <https://www.hdfgroup.org/solutions/hdf5>`_, a `matlab 
file <https://www.mathworks.com/products/matlab.html>`_, a plain text file but 
in a different format than what is required by 
:class:`~data_slicer.dataloading.Dataloader_3dtxt` or anything else.
In any case, you will have to find a way of loading that dataset into 
python first
(For the examples listed above, the packages `igor 
<https://pypi.org/project/igor/>`_, `h5py <https://github.com/h5py/h5py>`_, 
`scipy 
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.loadmat.html>`_ 
or `numpy 
<https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html>`_ could 
be used, respectively).

In order to make the example more concrete and to allow for a step-by-step 
code-along experience, we will explain here how we prepared the MRI brain scan 
data that can be loaded by issuing ``mw.brain()`` on ``PIT``'s console.
But again, how exactly this step is handled depends on your starting point.
The result of Step 1 should always be the same, though: Your data should be 
accessible as a numpy array in python.

OK, now let's get started with the concrete example.
We find and download our brain scan data of a meditating person 
`here <https://openneuro.org/crn/datasets/ds000108/snapshots/00002/files/sub-01:anat:sub-01_T1w.nii.gz)>`_ [#]_.

Once downloaded, we need a way of opening it.
A quick `Ecoisa search <www.ecosia.org>`_ leads us to the library `NiBabel 
<https://nipy.org/nibabel/gettingstarted.html>`_, which seems to be able to 
open ``.nii`` files.
Thus, we install that library (depending on your system and setup there may 
be different ways of doing this)::

   pip install nibabel

Now, following the instructions on the NiBabel webpage, we load the image 
data as a numpy array::

   python
   >>> import nibabel as nib
   >>> img = nib.load('sub-01_T1w.nii.gz')
   >>> my_data = img.get_fdata()

.. note::
   You need to be in the directory where you placed the downloaded file (here 
   ``sub-01_T1w.nii.gz``) in order for this to work.

.. note::
   Just to point it out once more, the details of this first step depend very 
   much on your use case. It also does not matter whether you work in the 
   live python interpreter like in the example or whether you wrap it all in 
   a script.
   You're fine as long as you have a way of getting your data into the form 
   of a numpy array.

Step 2: Optionally change data arrangement
..........................................

Now that we have our data in a numpy array, we are free to swap axes, cut off 
undesired parts or apply any processing we like to it.
This, again, depends completely on your use case.
Since in the example here we just want to *see* the data, we have nothing to 
do.  

In case you find yourself wanting to do some rearrangements, here are a few 
functions that might be of interest: :func:`numpy.moveaxis`, 
:func:`numpy.transpose` and 
all basic :class:`numpy.ndarray` operations, like `slicing and indexing 
<https://numpy.org/devdocs/reference/arrays.indexing.html#arrays-indexing>`_.

Step 3: Optionally create axis information
..........................................

Skipping this step means that the data we end up loading into ``PIT`` will 
have axes that simply count the number of pixels (voxels) from 0 upwards.
But we can assign more meaningful units to our axes, like in our 
example we could assign length units to the *x*, *y* and *z* axes.
To do this, we have to create one 1D array for each axis and collect them in 
a list.

In our example, we found by inspecting the original data that 1 pixel 
corresponds to 0.85 mm along the x and y directions and 1.5 mm along the z 
direction.
To create some reasonable axes, we could therefore do the following::

   >>> import numpy as np
   >>> nx, ny, nz = my_data.shape
   >>> x_axis = np.arange(0, nx*0.85, 0.85)
   >>> y_axis = np.arange(0, ny*0.85, 0.85)
   >>> z_axis = np.arange(0, nz*1.5, 1.5)
   >>> my_axes = [x_axis, y_axis, z_axis]

The three axes should of course have the lengths corresponding to the data 
dimensions.

.. 
   Usually your data will be a function of some variables, e.g. our brain 
   scan data is intensity as function of space coordinate :math:`I(x, y, z)`.
   Other things would be imaginable, for example pressure in the `xy` plane as a 
   function of time `t` :math:`p(x, y, t)` or intensity as a function of 
   momentum and energy :math:`I(k_x, k_y, E)`, etc.
   In our example with the brain, we are in the first mentioned situation 
   (:math:`I(x, y, z)`).

Step 4: Pickle it!
..................

Finally we can store our data in a format that can be efficiently read by 
``PIT``.
Here, we have different options, depending on whether or not we want to 
provide axes information (step 3).
In all three cases we make use of the convenience function 
:func:`~data_slicer.dataloading.dump`, which uses the ``pickle`` module to 
store any python object::

   >>> from data_slicer.dataloading import dump

Option 1: no axes information
*****************************

This is the easiest, you can just do::

   >>> dump(my_data, 'brain.p')

This will create the file ``brain.p`` in your current working directory.
If a file of that name already exists, it will ask you for confirmation.
(Obviously you can pick a filename of your choice. It doesn't even have to 
end in ``.p``.)

Option 2: with axes information in a dictionary
***********************************************

In order to also store the axis information we created in step 3, we just 
construct a :class:`dictionary <dict>` and pickle it::

   >>> D = dict(data=my_data, axes=my_axes)
   >>> dump(D, 'brain.p')

In this case it is important that the argument names ``data`` and ``axes`` 
are exactly like that. Other names will not work.
As the only exception, an alternative method is possible if you provide the 
three axes separately, like this::

   >>> D = dict(data=my_data, xaxis=x_axis, yaxis=y_axis, zaxis=z_axis)
   >>> dump(D, 'brain.p')

Option 3: with axes information in a Namespace
**********************************************

This option is given for convenience and out of consistency with the 
:class:`data_slicer.dataloading.Dataloader` objects.
Whether you use options 2 or 3 is entirely up to your personal preference and 
shouldn't make any difference.
The idea is exactly the same, except that we create a 
:class:`argparse.Namespace` instead of a dictionary::

   >>> from argparse import Namespace
   >>> D = Namespace(data=my_data, axes=my_axes)
   >>> dump(D, 'brain.p')

Conclusion
..........

And that's it. We have now successfully converted a datafile into a 
``PIT``-readable format.
Of course, if you have to do this kind of operation often, it would be a good 
idea to write a little script that does these steps for you.
If you're feeling confident, you could even create a :ref:`plugin 
<sec-plugin>` for the filetype(s) you need to use and make it available to 
other people.
Or, if you're lucky, somebody else has already done this and you can just use 
that plugin.

Plain text files
''''''''''''''''

Working with plain text (ASCII) files is significantly slower and requires 
more disk space than other file formats, but it can be useful to have the 
data in a human-readable form.
In order to create an plain text file in the correct format from some 
existing data, you will have to go through steps 1 to 3 exactly as in the 
description above.
The only thing that changes is the final step, step 4.

Step 4 for plain text files
...........................

In this case, we can just use the function 
:func:`data_slicer.dataloading.three_d_to_txt`::

   >>> from data_slicer.dataloading import three_d_to_txt
   >>> three_d_to_txt('brain.txt', my_data, axes=my_axes)

If you've skipped step 3, you can just leave out the ``axes`` argument.
In case you're typing along this tutorial, you will notice that the creation 
of this ``txt`` takes much longer than in the binary case - up to several 
minutes even.

.. rubric:: Footnotes

.. [#] This data set is taken from the OpenNeuro database.
       Openneuro Accession Number: ds000108
       Authored by: Wager, T.D., Davidson, M.L., Hughes, B.L., Lindquist, 
       M.A., Ochsner, K.N. (2008). Prefrontal-subcortical pathways mediating 
       successful emotion regulation. Neuron, 59(6):1037-50. 
       doi: ``10.1016/j.neuron.2008.09.006``

Examples
========

A couple of simple examples can be run to showcase some of the widgets in the 
package.
These are found in the ``tests`` submodule and you can try them out using::

   python -m data-slicer.tests.XXX

where ``XXX`` is the name of the test to run (name of the file in the 
``tests`` submodule, without ``.py`` ending):

   - ``threedwidget``:
     A simple widget with xy, xz and yz planes that can be moved through the 
     data cube.

   - ``freeslice``:
     shows a widget where the data cube can be sliced along the xy plane but 
     also freely with a Cutline.

   - ``pit``:
     just starts PIT.

Check out the source files at ``/your/package/location/data_slicer/tests/`` 
or directly on the `github 
<https://github.com/kuadrat/data_slicer/tree/master/data_slicer/tests>`_ to 
see how these widgets can be used.

.. _sec-access:

Accessing presented data
========================

Here follows a short overview on how to access some of the more useful data 
items represented in the different plots.
Please refer to the :ref:`graphical overview <sec-quickstart>` to understand 
which plot is being talked about.

.. note::
   Notice that accessing data in this way will *not* allow you to manipulate 
   the currently displayed data, i.e. you get a copy of whatever is currently 
   displayed and changes to this copy will not be reflected visually.  
   Only changes to the core data object (which is accessible through 
   ``pit.get_data()`` and settable through ``pit.set_data(...)``) and usage 
   of respective builtin functions (and likely functions from plugins) will be 
   visible.


Main and Cut Plot
-----------------

The currently displayed data of the *main plot* and *cut plot* can be 
accessed respectively by::

   pit.get_main_data()
   pit.get_cut_data()

or alternatively::

   mw.main_plot.image_data
   mw.cut_plot.image_data

Profiles
--------

The *horizontal* and *vertical profiles* as well as the *integrated 
intensity* profile can be accessed by::

   pit.get_hprofile()
   pit.get_vprofile()
   pit.get_iprofile()

.
.. _sec-console:

Using the console
=================

An essential part of PIT is the console.
This is just a regular `ipython console <https://ipython.org/>`_ which allows 
you, in principle, to do anything you could imagine doing with python:
set and use variables, define and run functions, run scripts, use matplotlib 
to create plots, etc.
Check out the :ref:`ipython crash course <sec-ipython-crash-course>` if you 
are unfamiliar with ipython (or python).

Importantly, all the visible objects and the underlying data structures can 
be accessed through the console.
There are two objects through which this access is provided: ``pit`` and ``mw``.

``mw`` is short for *main window* and contains all the visible elements, e.g. 
the different plots (``mw.main_plot``, ``mw.cut_plot``, ``mw.x_plot``, 
``mw.y_plot``, etc.)
The colormap setting is also handled through ``mw``, through 
:func:`mw.set_cmap() <data_slicer.pit.MainWindow.set_cmap>`.

``pit`` is an instance of :class:`~data_slicer.pit.PITDataHandler` and is 
responsible for keeping all the visible and invisible data elements consistent.
Use it to load data in (:func:`pit.open() <data_slicer.pit.PITDataHandler>`), 
change the orientation in which we look at the data (:func:`pit.roll_axes() 
<data_slicer.pit.PITDataHandler.roll_axes`).
Here's an incomplete list of some of the more useful functions of the ``pit`` 
object:

   - **data access**

     Confer :ref:`sec-access` for details on how to access various data 
     elements.

   - :func:`pit.reset_data() <data_slicer.pit.PITDataHandler.reset_data>`

      Reset everything back to the state just after loading the data in.

   - :func:`pit.lineplot() <data_slicer.pit.PITDataHandler.lineplot>`

     Create a matplotlib figure that shows the data as a series of lines.

   - :func:`pit.plot_all_slices() 
     <data_slicer.pit.PITDataHandler.plot_all_slices>`

     Show equidistantly spaced slices along z in matplotlib figures.

   - :func:`pit.overlay_model() <data_slicer.pit.PITDataHandler.overlay_model>`
    
     Supply a function of the two variables (x and y axes of the *main plot* 
     and display it in the *main* and *cut plots*.

   - :func:`pit.remove_model() <data_slicer.pit.PITDataHandler.remove_model>`
     
     Remove the lines from above command.

.. _sec-quickstart:

Quick start
===========

If you want to dive right in, just type ``pit`` from a command line.
The startup data is loaded and you can familiarize yourself with the layout 
and the most basic functionality.

You can alternatively load a set of MRI brain scan data that is perhaps more 
intuitive to understand by using the command ``mw.brain()`` on the console 
[#]_.

.. figure:: ../../img/pit_overview.png
   :scale: 50 %
   :alt: Image not found.

   The main window of PIT.
   
   =  ==========================================================================
   a  main data plot; ``mw.main_plot``
   b  cut plot; ``mw.cut_plot``
   c  vertical profile; ``mw.x_plot``
   d  integrated z plot; ``mw.integrated_plot``
   e  horizontal profile; ``mw.y_plot``
   f  colorscale sliders
   g  interactive ipython console
   =  ==========================================================================

Main data plot and slice selection
----------------------------------

If you imagine your 3D dataset as a cube (*data cube*), the *main data plot* 
**(a)** would initially represent what you would see when looking at the cube 
from the top.
The horizontal **x** axis corresponds to the first dimension of your data 
cube, while the vertical **y** axis corresponds to the second.
The *integrated z plot* **(d)** shows the sum of each xy slice along the 
third **z** dimension.
You can drag the yellow slider to select a different slice to be displayed in 
the *main data plot*.
Additionally, using the **up** and **down arrow keys** (after having clicked in 
the *integrated z plot*) allows you to in- or decrease the integration range 
for the slice.
**Left** and **right arrow keys** move the slider step by step.


Creating arbitrary cuts
-----------------------

The yellow line inside the *main data plot* is called the *cutline*.
It is draggable as a whole and at the handles.
This is the knife that cuts through our data cube and the *cut plot* **(b)** 
shows what we see when we cut the data cube along the cutline.
Hitting the **r key** will re-initialize the cutline alternatingly at an 
angle of 0 or 90 degrees.
This is also useful if you happen to "lose" the cutline.

The *horizontal* and *vertical profiles* **(c)** and **(e)** just display the 
line profiles of the data shown in the *cut plot* along the cursor.


Colorscale sliders and the ipython console
------------------------------------------

The *colorscale sliders* **(f)** enable you to quickly change the min and max 
values of the colorscale as well as the exponent of the powerlaw 
normalization (*gamma*).

To change to used colormap, the *ipython console* **(g)** has to be used:
Type the command ``mw.set_cmap('CMAP_NAME')``.
Refer to :ref:`the section about using the console <sec-console>` for more.


Basics
------

Data is loaded by using the :func:`pit.open() 
<data_slicer.pit.PITDataHandler.open>` command from the console.
Change the way we look at the data using :func:`pit.roll_axes() 
<data_slicer.pit.PITDataHandler.roll_axes>` and overlay a model over the 
displayed data with :func:`pit.overlay_model() 
<data_slicer.pit.PITDataHandler.overlay_model>`.

You can create a matplotlib figure of the *main* or *cut plots* by right 
clicking and choosing *MPL Export*.

.. note::
   The *Export...* option from the right click menu is broken due to some 
   error in pyqtgraph over which I have no control.

Refer to :ref:`sec-console` for more information on what you can do.

.. rubric:: Footnotes

.. [#] This data set is taken from the OpenNeuro database.
       Openneuro Accession Number: ds000108
       Authored by: Wager, T.D., Davidson, M.L., Hughes, B.L., Lindquist, 
       M.A., Ochsner, K.N. (2008). Prefrontal-subcortical pathways mediating 
       successful emotion regulation. Neuron, 59(6):1037-50. 
       doi: ``10.1016/j.neuron.2008.09.006``

.. _sec-config:

The config directory
====================

You can customize a few things about PIT on your system by editing its
configuration directory.
You can find the location of this directory by issuing::

   pit.get_config_dir()

Colormaps
---------

In the ``cmaps`` subdirectory under your config dierectory (you can create it 
if it does not exist) you can put custom colormaps that will then be 
available to PIT.
Colormaps should be simple text files with three columns, representing RGB 
values (the scale doesn't matter), so it could look like this::

    0 0 0
    200 0 0
    400 0 0
    800 0 10
    [...]

The filename (everything before the first ``.`` that appears in the filename) 
will be the name under which you can find and load that colormap from within 
PIT.

.. _sec-config-plugins:

Plugins
-------

Any python module that you put in the ``plugins/`` subdirectory will be found 
by PIT even if it is not part of the system PYTHONPATH.

Additionally, if you create a simple txt file called ``autoload.txt`` in 
which you list plugins (the names you would use to import them, one per line), 
PIT will attempt to load all of them on startup.

Example contents of an ``autoload.txt`` with just one plugin::

   ds_arpes_plugin

.. _sec-ipython-crash-course:

ipython crash course
====================

For a full tutorial on ipython, visit its `offcial documentation 
<https://ipython.readthedocs.io/en/stable/interactive/index.html>`_.
The pages here are just to point out the most useful features.

Basically, ipython is an improved python interpreter.
Just like in the usual python interpreter, you can execute any python code 
you want.
Code you enter can even span multiple lines, allowing you to define functions 
or classes::

   In [1]: print('Hello World.')
   Hello World.

   In [2]: def my_function(a) :
      ...:     return a**2
      ...:     

   In [3]: my_function(2)
   Out[3]: 4

One of the most useful features of ipython is the tab-completion: pressing 
Tab while writing some code in ipython will make the interpreter try to guess 
what you want to type and complete it.
If it isn't sure, it will present you with a list of options that can be 
scrolled through using Tab or the arrow keys.
Try the following as an example (everything after the ``#`` should not be 
typed, but understood as a comment)::

   In [1]: my_long_variable = 3.14

   In [2]: my # press Tab now and it will complete to `my_long_variable`
   Out[2]: 3.14

   In [3]: my_other_long_variable = 2.72
   
   In [4]: my # press Tab now and you will be presented with a list 
              # containing `my_long_variable` and `my_other_long_variable`

This is not just convenient to save some typing, but even more importantly, 
it helps you check what options you have.
In PIT, you might forget the names of different functions.
In such a situation you could use Tab to get a list of all functions on the 
fly::

   In [1]: pit. # press Tab and you will get a list of all functions the 
                # ``pit`` object offers

You are able to do that with any python object, and even with filenames (e.g. 
when you're doing a ``pit.open('...')``).

Another useful features are the ipython *magic* commands::

   In [1]: %run my_script.py # This runs the file `my_script.py`, if it is in 
                             # the current working directory. All the 
                             # variables and objects created in that script 
                             # will afterwards be available to you in the 
                             # ipython session
   Out[1]: ...

   In [2]: %load my_script.py # This will load the contents of `my_script.py` 
                              # into the interpreter, as if you had typed it 
                              # in by hand. This would allow you to review 
                              # and makes changes to the code before running them.

Do not forget the existence of python's built in ``help()`` function, to 
quickly get reminded of how a function should be used without having to search online.

This should sum up the most useful features for the purposes of using PIT.
Again, for more details on ipython, visit `the offcial documentation 
<https://ipython.readthedocs.io/en/stable/interactive/index.html>`_.

.. _sec-plugin:

Plugins
=======

It is possible for people to write plugins that add additional functionality 
to PIT (see also :ref:`the plugin writing guide <sec-plugin-how-to>`).
A plugin for PIT is just a python module that connects to the backside of PIT.
This module should either be found by the python version you are running 
(i.e. you should be able to ``import`` it) or you can put it in your 
data-slicer :ref:`config directory <sec-config>`.

Loading plugins
---------------

If you have your plugin installed or ready in your config directory, you can 
load it from PIT by issuing::
   
   my_plugin = mw.load_plugin('PLUGIN_NAME')
   
Of course you can use any variable name in place of the example's 
``my_plugin``. ``PLUGIN_NAME`` should be the exact same thing you would 
write when you import the module using python's ``import`` statement.  

Autoloading plugins
-------------------

You can set it up such that your favourite plugins will always be loaded on 
startup of PIT.
See :ref:`here <sec-config-plugins>`.

.. _sec-plugin-list:

List of known plugins
---------------------

If you have created a PIT plugin, feel free to add it to the list below, 
either via pull request or through an e-mail 
(kevin.pasqual.kramer@protonmail.ch).

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Plugin Name (link) 
     - Description

   * - `ds-example-plugin <https://github.com/kuadrat/ds_example_plugin>`_
     - exists as a minimal example that can be used for guidance when 
       creating your own plugins. A step-by-step tutorial on how it was made 
       can be found :ref:`here <sec-plugin-how-to>`.

   * - `ds-arpes-plugin <https://github.com/kuadrat/ds_arpes_plugin>`_
     - tools for angle-resolved photoelectron spectroscopy (ARPES) data

Installation
============

The installation of ``data-slicer`` has been tested on Linux, macOS and Windows.

The easiest way to install the package is to use 
`pip <https://pip.pypa.io/en/stable/>`_. Just type the following on a command 
line::

   pip install data_slicer

If the above does not appear to work, it might help to try installation from 
a virtual environment. 

Anaconda
--------

Detailed instructions for Anaconda users follow:

1) Open "Anaconda Prompt" 

2) In order not to mess up your dependencies, create a virtual 
   environment with python version 3.7.5::

      $ conda create --name testenv python==3.7.5
      [some output]
      $ conda activate testenv
      (testenv) $

3) Inside your virtual environment, run the following commands to download and 
   install data_slicer with all its dependencies (the first one is just to 
   upgrade pip to the latest version)::
   
      (testenv) $ pip install --upgrade pip
      (testenv) $ pip install data_slicer
   
   This will create a lot of console output. If everything succeeded, you should 
   see something like ``Successfully installed data_slicer`` towards the end.

4) Test the installation by starting up PIT::

      (testenv) $ pit
   
   This should bring up a window with some example data.
   See also :ref:`below, for how to run automated tests to verify your 
   installation <sec-verifying>`.

.. _sec-verifying:

Verifying your installation
---------------------------

Once installed, you can run a set of automated tests in order to check if the 
main features work as expected.
To do this, issue the following on the command line::

   python -m data_slicer.tests.run_tests

The result should be that some text is printed to the console and some 
windows open, with a few things happening in them before they quickly close 
again.
Basically, these tests simulate a few interactions that the user could have 
with these windows and verify that they worked with some checks.
If all went well you might see some warnings, but no notifications of any
``failures``.
It could, for example, look like this::
   
   ================== 4 passed, 16 warnings in 14.92 s ==================

.. note::
   The fact that all tests passed does not guarantee that *everything* is in 
   working order - but it's a very good sign.

If interested, you can also run these tests individually and interact with 
the respective windows by calling them like so::

   python -m data_slicer.tests.test_XXX

where ``XXX`` is any of ``pit``, ``freeslice``, ``threedwidget``.

Upgrading
---------

The following command will attempt to upgrade ``data-slicer`` to the latest 
published version::

   pip install --upgrade data_slicer

It is usually a good idea to upgrade ``pip`` itself before running above 
command::

   pip install --upgrade pip

.. Note::
   Run these commands from within the same (virtual) environment as you've 
   installed ``data-slicer`` in.

Dependencies
------------

This software is built upon on a number of other open-source frameworks.
The complete list of packages is:

.. include:: ../../requirements.txt
   :code:

Most notably, this includes 
`pyqtgraph <https://pyqtgraph.readthedocs.io/en/latest/>`_ for fast live 
visualizations and widgets, `numpy <https://numpy.org/>`_ for numerical 
operations and `matplotlib <https://matplotlib.org/>`_ for plot exporting 
functionalities.

.. _sec-widgets:

Widgets
=======

The widgets shipped with ``data-slicer`` can be re-combined to create new 
applications.
In order to do this, a certain level of familiarity with widget-based GUI 
creation is required and some knowledge or experience with `Qt 
<https://www.qt.io/>`_ or `PyQt <https://wiki.python.org/moin/PyQt>`_ is 
recommended.
Nevertheless, even without said experience, the `step-by-step 
<sec-widget-tutorial>`_ tutorial below may help you figure out how to do such 
things.

List of available widgets
-------------------------

.. list-table:: 
   :header-rows: 0
   :widths: 40 60

   * - :class:`data_slicer.imageplot.ImagePlot`
     - A pseudocolor plot of a 2D dataset, similar to :class:`matplotlib 
       pcolormesh <matplotlib.pyplot.pcolormesh>`.

   * - :class:`data_slicer.imageplot.CrosshairImagePlot`
     - An :class:`~data_slicer.imageplot.ImagePlot` with a draggable crosshair.

   * - :class:`data_slicer.imageplot.CursorPlot`
     - A regular data plot with a draggable line (cursor).

   * - :class:`data_slicer.imageplot.Scalebar`
     - A draggable scalebar. 

   * - :class:`data_slicer.cutline.Cutline`
     - A line that can be dragged on both ends and added to a
       :class:`pyqtgraph.PlotWidget` to create arbitrary cuts.  

   * - :class:`data_slicer.pit.MainWindow`
     - The full main window of PIT, itself consisting of widgets from this 
       list.  

   * - :class:`data_slicer.widgets.ColorSliders`
     - Gamma and vmax color sliders. Essentially just a pair of 
       :class:`~data_slicer.imageplot.Scalebar` objects.

   * - :class:`data_slicer.widgets.ThreeDWidget`
     - A widget that allows showing xy-slices out of a 3D dataset as a plane in 
       3D.

   * - :class:`data_slicer.widgets.ThreeDSliceWidget`
     - Subclass of :class:`~data_slicer.widgets.ThreeDWidget` that shows the 
       xz- and yz-planes in addition to the xy plane.

   * - :class:`data_slicer.widgets.FreeSliceWidget`
     - Subclass of :class:`~data_slicer.widgets.ThreeDWidget` that makes use 
       of a :class:`~data_slicer.cutline.Cutline` on an 
       :class:`~data_slicer.imageplot.ImagePlot` to allow creating arbitrary 
       slice-planes.

.. _sec-widget-tutorial:

Example
-------

In the following we will go through the creation of a simple App that uses a 
:class:`~data_slicer.widgets.ThreeDSliceWidget` to display some data that can 
be loaded when clicking a **Load** button.

Step 1: Just a plain ThreeDSliceWidget
......................................

First, we take a look at what we would need to do in order to just create a 
:class:`~data_slicer.widgets.ThreeDSliceWidget` without anything else.

.. literalinclude:: ../examples/example_step1.py
   :linenos:

This should create a window with a 
:class:`~data_slicer.widgets.ThreeDSliceWidget` - but there is no data 
present inside the widget yet.

Step 2: Adding a Load button
............................

In order to add a button, we need to change the structure of our application 
a little bit.
Since we need to be able to let Qt know where our different GUI elements 
should appear, we are going to work with a :class:`~Qt.QtGui.QGridLayout`.
We then add a :class:`~Qt.QtGui.QPushButton` and the 
:class:`~data_slcier.widgets.ThreeDSliceWidget` from before to the layout.
In the following code snippets, all new lines are preceded by a comment ``# 
NEW``.

.. literalinclude:: ../examples/example_step2.py
   :linenos:

We now have a button above the :class:`~data_slcier.widgets.ThreeDSliceWidget`!
However, clicking the button does not do anything yet... Let's change that in 
the next step.

Step 3: Making the button do something
......................................

We can define what happens when the button is clicked by *connecting* a 
function to it::

   load_button.clicked.connect(my_function)

For ``my_function`` we could put any ``python`` callable, for example::

   def my_function() :
      print('Button clicked!;)

Try defining ``my_function`` like this **before** connecting it to the button 
and running the example again.
You should now get the message ``Button clicked!`` on the console  whenever 
you click the button.

This is not yet what we want though.
We would like the click on the **Load** button to open a file selection 
dialog from which we can navigate to a data file, select it and that is then 
going to be loaded into the :class:`~data_slicer.widgets.ThreeDSliceWidget`.
This can be achieved with the following function::

   def load_data() :
       # Open the file selection dialog and store the selected file as *fname*
       fname = QtGui.QFileDialog.getOpenFileName(layout_widget, 'Select file')
       print(fname[0])
       
       # Load the selected data into the ThreeDSliceWidget
       D = dataloading.load_data(fname[0])
       widget.set_data(D.data)

Don't forget to ``from data_slicer import dataloading`` at the start of the 
file and to connect this function to our **Load** button.

The full example code at this stage should look like this:

.. literalinclude:: ../examples/example_step3.py
   :linenos:

Conclusion
..........

While we have seen how to use the provided widgets in other contexts to 
create new applications, it is obvious that a lot could be improved and 
tinkered with in our little example.
One could include some :class:`~data_slicer.widgets.ColorSliders` and link 
them up to the :class:`~data_slicer.widgets.ThreeDSliceWidget` s colormap to 
adjust the colors.
Or one could implement a different way of loading the data that is not 
limited to the formats supported by :mod:`data_slicer.dataloading`, add 
error handling when specifying unsupported formats, add more widgets that 
show the data from different viewpoints, and so much more (don't even get me 
started on brushing up the layout and *look and feel*).

For further inspiration it is recommended to check out the source coudes of 
the tests, located under the ``tests`` directory in your ``data_slicer`` 
installation or on `the github 
<https://github.com/kuadrat/data-slicer/tree/master/data_slicer/tests>`_.
Additionally, a lot can be achieved when combining functionalities from 
`pyqtgraph <https://pyqtgraph.readthedocs.io/en/latest/>`_.
Check out the rich set of examples they provide by::

   python -m pyqtgraph.examples
Welcome to data-slicer's documentation!
=======================================

Welcome to the documentation for the ``data-slicer`` package!
These pages are intended to give an overview about what the package provides, 
and explain how it can be installed and used.

.. toctree::
   :maxdepth: 1
   :caption: Overview

   introduction
   installation

.. toctree::
   :maxdepth: 1
   :caption: PIT

   quickstart
   data_format
   console
   access

.. toctree::
   :maxdepth: 1
   :caption: Other

   plugins
   configuration
   ipython_crash_course
   widgets
   examples
   contributing

.. toctree::
   :maxdepth: 1
   :caption: Code reference

   data_slicer

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Introduction
============

The ``data-slicer`` package provides tools for quick visualization of 3 
dimensional datasets.
These tools are not intended to enable the creation of publication quality 
renderings, but rather to allow users to quickly inspect datasets in order to 
get an overview of the data and to quickly visualize the result of analytical 
operations on the data under study.
Examples where this package has been and could be employed included, but are 
not limited to:

   - medical or biological imaging techniques.

   - synchrotron based scattering and spectroscopic experiments, where 
     deciding on the course of the time-limited experiment often requires 
     preliminary analysis.

   - visualizing satellite imagery or georeferenced data.

Basically, ``data-slicer`` can be employed in all scientific, engineering, 
medical, artistic or other data driven disciplines where the use case of 

inspecting and slicing data (hyper)cubes arises.

At its core, the package contains the Python Image Tool, PIT for short, a 
simple but fast GUI to efficiently look at 3D data, slice it in various ways 
and employ data processings from the command line of an ``ipython`` console. 
While being set up to work in a general manner, specific processing and 
analysis routines for different fields can be integrated into PIT by means of 
a simple :ref:`plugin mechanism <sec-plugin>`.

Besides PIT, other widgets for 3D visualization are also available.
Thus, the package is intended to provide user bases in different fields with 
a means of easily creating customized tools for their data inspection and 
analysis routines.

data\_slicer package
====================

The following is an auto-generated overview over the full package structure.
You should be able to find some help about every module, class and function 
contained in the `data_slicer` package.

PIT and MainWindow (the data\_slicer.pit module)
------------------------------------------------

These classes are what you will be mostly interacting with when using PIT.
Of most interest should be :class:`data_slicer.pit.PITDataHandler` and 
:class:`data_slicer.pit.MainWindow` as these two are accessible from PIT's 
ipython console as ``pit`` and ``mw``.

.. automodule:: data_slicer.pit
   :members:
   :undoc-members:
   :show-inheritance:

Other modules
-------------

The following modules constitute the base of `pit` and other tools provided 
by `data_slicer`.

data\_slicer.cmaps module
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.cmaps
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.cutline module
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.cutline
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.dataloading module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.dataloading
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.dsviewbox module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.dsviewbox
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.imageplot module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.imageplot
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.model module
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.model
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.plugin module
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.plugin
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.set\_up\_logging module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.set_up_logging
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.utilities module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.utilities
   :members:
   :undoc-members:
   :show-inheritance:

data\_slicer.widgets module
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: data_slicer.widgets
   :members:
   :undoc-members:
   :show-inheritance:

