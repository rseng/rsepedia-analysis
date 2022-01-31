# Pylustrator

[![DOC](https://readthedocs.org/projects/pylustrator/badge/)](https://pylustrator.readthedocs.io)
[![Build Status](https://app.travis-ci.com/rgerum/pylustrator.svg?branch=master)](https://app.travis-ci.com/github/rgerum/pylustrator)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0.html)
[![DOI](https://img.shields.io/badge/DOI-10.21105/joss.01989-blue.svg)](https://doi.org/10.21105/joss.01989)

<img style="float: left;" alt="docs/images/logo.png" src="docs/images/logo.png" />

Pylustrator is a software to prepare your figures for publication in a reproducible way. This means you receive a figure
representing your data and alongside a generated code file that can exactly reproduce the figure as you put them in the
publication, without the need to readjust things in external programs.

Pylustrator offers an interactive interface to find the best way to present your data in a figure for publication.
Added formatting an styling can be saved by automatically generated code. To compose multiple figures to panels,
pylustrator can compose different subfigures to a single figure.

Please also refer to the [Documentation](https://pylustrator.readthedocs.io) for more information.

## Issues, Questions, and Suggestions

Please submit your questions, suggestions, and bug reports to the
[Issue Tracker](https://github.com/rgerum/pylustrator/issues)


## Contributing

You want to contribute? Great!
Contributing works best if you creat a pull request with your changes.

1. Fork the project.
2. Create a branch for your feature: `git checkout -b cool-new-feature`
3. Commit your changes: `git commit -am 'My new feature'`
4. Push to the branch: `git push origin cool-new-feature`
5. Submit a pull request!

If you are unfamilar with pull requests, you find more information on pull requests in the
 [github help](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests)
---
title: 'pylustrator: code generation for reproducible figures for publication'
tags:
  - reproducibility
  - code generation
  - interactive
  - Python
  - matplotlib
  - plotting
  - drag
  - style
authors:
  - name: Richard Gerum
    orcid: 0000-0001-5893-2650
    affiliation: 1
affiliations:
 - name: Department of Physics, University of Erlangen-Nürnberg, Germany
   index: 1
date: 06 September 2019
bibliography: paper.bib
---


# Background

In recent years, more and more researchers have called attention to a growing "reproducibility crisis" [@CRL16846]. An important factor contributing to problems in the reproducibility of results from published studies is the unavailability of the raw data from the original experiment and the unavailability of the methods or code used for the analysis of raw data [@Baker2016]. One major step to overcome these shortcomings is publishing all raw data and a documented version of the code used for analysis [@baker2016scientists]. Ideally, anyone interested should be able to download the raw data and reproduce the figures of the publication exactly.

To address the issue of data availability, researchers are encouraged to make their data available in online repositories such as Dryad [@dryad]. However, these data are useless unless the complete analysis procedure, including all analysis and visualisation steps, can be comprehended by other scientists. The best way to achieve this is to provide a complete, well-documented analysis code, including all important steps from the basic artifact corrections to the final plot to be published. Open source scripting languages such as Python [@Rossum1995] or R [@R] are ideal for such code because open source languages are accessible to everyone. In addition, interpreted languages do not need to be compiled, therefore present fewer obstacles for the user to run the code. The final part of the data analysis is the visualisation, which is crucial for communicating the results [@Tufte1893]. This paper deals with the visualization step, consisting of two parts: generating simple plots from data, and composing meaningful figures from these plots.

The first part of generating the building blocks of figures, the plots, is already covered in various toolkits, e.g., Matplotlib [@Hunter2007], Bokeh [@Bokeh] or Seaborn [@seaborn]. There already exist some software that records user interactions to generate 3D plots [@ahrens2005paraview,@westenberger2008avizo,@dragonfly], but no convenient Python toolkit is yet available to generate reproducible figures from simple plots scripts. Matplotlib offers figures composed of several subplots, but to create a complete, publication-ready figure a lot of code is needed to add all formatting, annotation and styling commands. Users often prefer graphical tools like image manipulation software, e.g., GIMP [@GIMP] or Inkscape [@Inkscape]. These offer great flexibility, but do not provide a reproducible way of generating figures and bear the risk of accidentally changing data points. It is also important to note that when using an image manipulation software, every small change in the analysis requires re-editing the figure in the image manipulation software. This process slows down the creation of figures and is prone to errors.
The ``pylustrator`` package was developed to address this issue.

# Algorithm and Examples

``Pylustrator`` fills the gap from single plots to complete figures via a code-generation algorithm that converts user input into Python code for reproducible figure assembly (Fig. 1).  Minor changes to the analysis or new data only require re-execution of the code to update the figure.

![Example of composing a figure with pylustrator.](figure1.pdf)

Using ``pylustrator`` in any Python file that uses Matplotlib to plot data requires only the addition of two lines of code:

    import pylustrator
    pylustrator.start()
    
The Matplotlib figure is then displayed in an interactive window (Fig. 2) when the command `plt.show()` is called. In this interactive window ``pylustrator`` allows the user to:

- resize and position plots by dragging with the mouse 
- adjust the position of plots legends
- align elements easily through automatic "snap-in" mechanism
- change the size of the whole figure in cm/inch
- add text and annotations and change their style and color
- adjust plot ticks and tick labels

![The interface of ``pylustrator``. The user can view the elements of the plot, edit their properties, edit them in the plot preview and experiment with different color schemes.](figure2.pdf)

``pylustrator`` tracks all changes to the figure and translates them into Python code. 
To do so, the internal representation of changes has to fulfill some requirements. 
Changes need to be able to be replaced by newer changes that affect the same property of the same object, they
need to be able to be converted to code, and changes need to be retrieved from generated code when loading a file that has already ``pylustrator``-generated code in it.

Each change is defined by two parts: the affected object (e.g., a text object) and the affected property 
(e.g., its color). If a change has the same object and property as a previous change, it overwrites the previous change.

Changes are converted to code by first serializing the affected object by iteratively going up the parent-child tree
 from, e.g., a text object, to the axis that contains the text to the figure that contains the axis. From this dependency relation, a Python code segment is generated (e.g., `plt.figure(1).axes[0].texts[0]`, the first text of the first axis of figure 1).
  Then the property command is added (e.g., `.set_color("#ff0000ff")`).
  When saving, ``pylustrator`` introspects its current execution stack to find the line of code from where it was called and inserts the automatically generated code directly before the command calling ``pylustrator``.
   
 When a file with automatically generated code is loaded (see code example in figure 1), ``pylustrator`` splits all the automatically generated lines into the affected objects and affected properties. New changes, where both the affected object and the affected property match a previous change, overwrite the previous change. This ensures that previously generated code can be loaded appropriately, and saving the same figure multiple times does not generate the line of code for this change multiple times.
  
It is important to note that the automatically generated code only relies on Matplotlib and does not need the ``pylustrator`` package anymore. Thus, the ``pylustrator`` import can later be removed to allow sharing the code without an additional introduced dependency. 

The documentation of ``pylustrator`` can be found on https://pylustrator.readthedocs.org.

# Conclusion
This packages offers an improvement to create publishable figures from single plots based on an open source Python library called ``pylustrator``. The figures can be arranged by drag-and-drop, and the ``pylustrator`` library generates the corresponding code. This library provides a valuable contribution to improve reproducibility of scientific results.

# Acknowledgements 

We acknowledge testing, support and feedback from Christoph Mark, Sebastian Richter, and Achim Schilling and Ronny Reimann for the design of the Pylustrator Logo.

# References
.. _composing:

Composing Figures
=================

Pylustrator can also be used to compose panels of different subplots using the function the function `pylustrator.load()`.

Example
-------

Suppose we have two plot files that we want to include in our figure:

.. sidebar:: Plot 1

    .. image:: plot1.png
       :height: 150px

.. literalinclude:: plot1.py
   :caption:
   :language: python
   :linenos:

.. raw:: html

    <div style="clear: both;"></div>

.. sidebar:: Plot 2

    .. image:: plot2.png
       :height: 150px

.. literalinclude:: plot2.py
   :caption:
   :language: python
   :linenos:

.. raw:: html

    <div style="clear: both;"></div>

Now we can create a script that generates the composite figure:

.. literalinclude:: figure1.py
   :caption:
   :language: python
   :linenos:
   :emphasize-lines: 4,5

.. image:: figure1.png

The subfigures can both contain also multiple axes (e.g. subplots) or be previously styled with pylustrator. The composite
figure can also be styled with pylustrator to finialize the figure.

.. Note:: Please note that the code of the target script will be executed, to only load script files from trusted sources.
    This also holds true when loading scripts from a cached pickle representation, see `Caching`_.

Supported Formats
-----------------

With `pylustrator.load()` different types of inputs can be used to add to a composite figure.

Python Files
~~~~~~~~~~~~
If the input is a ".py" file then the file is compiled and executed. All the figure elements that are generated from the
python file are then inserted into the target figure.

.. Note:: All your plotting code should be in the main part of the script, everything hidden behind a
    `if __name__ == "__main__":` will be ignored.

Image Files
~~~~~~~~~~~
If the input is any image format that can be opened with `plt.imread()` then the file is loaded and its content is displayed
in a separate axis using `plt.imshow()`. To scale the image, pylustrator uses the dpi settings of the figure. You can also
provide the `dpi` keyword to `pylustrator.load()` to scale the image.

Svg Files
~~~~~~~~~
If the input file is a `svg` file, then pylustrator tries to generate matplotlib patches to display the content of the svg
file in a new axes.

.. Warning:: The svg specification is broader than the patches supported by matplotlib. Therefore, gradients, filters,
    fill patterns and masks cannot be supported. Also svg offers some advances positioning features for text (e.g.
    letter-spacing, word-spacing) which are difficult to match in matplotlib. There might be some more differences in
    details in the implementations. If you thing something can be addressed in matplotlib, you can report it in the
    `bugtracker <https://github.com/rgerum/pylustrator/issues>`_.


Positioning
-----------
Loaded subfigures can be positioned using the `offset` keyword. Offsets can have different units.

The default unit is to
interpret it as a percentage of the current figure size. Note, that when calling multiple imports, the figure size changes
between calls. Therefore, a 2x2 grid can be achieved with relative offsets as follows:

.. code-block:: python
    :linenos:

    pylustrator.load("plot1.py")
    pylustrator.load("plot2.png", offset=[1, 0])
    pylustrator.load("plot3.jpg", offset=[0, 1])
    pylustrator.load("plot4.svg", offset=[0.5, 0.5])

The offset can also be specified in cm or in in. Therefore, the third entry of the tuple can be the string "cm" or "in".
For example:

.. code-block:: python
    :linenos:

    pylustrator.load("plot1.py")
    pylustrator.load("plot2.png", offset=[2, 1, "cm"])
    pylustrator.load("plot3.jpg", offset=[0.3, 0.9, "in"])

Caching
-------
Pylustrator also offers the possibility to cache the figures generated by a script file. Therefore, the figure is pickled
after it has been created from the script and saved next to the script. If the script is newer (last modified timestamp)
as the pickle file with the cached figure, the script is executed again. The caching behavior can be disabled with the
keyword `cache=False`.
.. _styling:

Styling Figures
===============

Opening
-------
To open the pylustrator editor to style a figure, you just have to call the function `pylustrator.start()` before any figure is
created in your code.

.. code-block:: python
    :linenos:

    import pylustrator
    pylustrator.start()

This will overload the commands `plt.figure()` and `plt.show()`. `plt.figure()` (which is often indirectly called when
creating plots) now initializes a figure in a GUI window and `plt.show()` then shows this GUI window.

In the GUI window elements can be dragged around using a click and drag with the left mouse button. If you want to cycle
though different elements on the same spot double click on the same position multiple times. To zoom in in the plot window
use `strg+mousewheel` and you can pan the figure with holding the middle mouse button.

To select multiple elements hold shift while clicking on multiple elements.

Saving
------

To save the figure press `ctrl+s` or select `File->Save`. This will generate code that corresponds to the changes you made
to the figure and past it into your script file or your jupyter notebook cell. The code will be pasted directly over the
`plt.show()` command that started the editor or, if there already is a generated code block for this figure, it will replace
the existing code block.

Increasing Performance
----------------------
Often plots with lots of elements can slow down the performance of pylustrator as with every edit, the whole plot is
rerendered. To circumvent this problem, pylustrator offers a function to calculate a rasterized representation (e.g. pixel data)
of the contents of each axes and only display the pixel data instead of rendering the vector data with every draw.

It can be activated with the button "rasterize". It can be clicked again to update the rasterisation or deactivated with
a click on the button "derasterize" next to it.

Color editor
------------
Pylustrator comes with a powerful color editor which allows to test different color configurations for your figure easily.
On the right hand side of the window you see a list of all currently used colors. You can right click on any color to open
a color choosed dialog. You can also directly edit the color using the html notation (e.g. #FF0000) provided on the button.
You can drag and drop colors to different slots to test different configurations for your figure.

The field below allows to copy and paste color lists from different sources. For example using color palette generators on the
internet, e.g. `<https://medialab.github.io/iwanthue/>`_.

Additionally, if you generate plot lines with colors from a colormap, pylustrator can recognize that and allow you to
choose different colormaps for the set of plot lines.

Tick Editor
-----------

.. |the tick icon| image:: ../pylustrator/icons/ticks.ico

To edit the ticks of an axes click on |the tick icon|. There, a windows opens that allows to set major and minor ticks
every line in the edit window corresponds to one tick. Texts are directly interpreted as float values if possible and the
text used as tick label, e.g. you can but a tick at 5.0 with the text "5" (e.g. formated without decimal point).
To specify an exponent, it is also possible to write e.g. 5*10^2 (to put a tick at 500 with the label :math:`5\cdot10^2`).
If the label cannot be directly writen as a number, add the label after the number enclosed in tick marks e.g. 5 "start",
to add a tick at position 5 with the label "start".

.. figure:: images/logo.png
    :align: left

Welcome to the Pylustrator Documentation
========================================

Pylustrator is a software to prepare your figures for publication in a reproducible way. This means you receive a figure
representing your data and alongside a generated code file that can exactly reproduce the figure as you put them in the
publication, without the need to readjust things in external programs.

Pylustrator offers an interactive interface to find the best way to present your data in a figure for publication.
Added formatting an styling can be saved by automatically generated code. To compose multiple figures to panels,
pylustrator can compose different subfigures to a single figure.


.. raw:: html

    <div style="clear:both"></div>
    <hr>

Installation
------------

Just get pylustrator over the pip installation:

    ``pip install pylustrator``

The package depends on:

numpy, matplotlib, pyqt5, qtpy, qtawesome, scikit-image, natsort

Usage
-----

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="//www.youtube.com/embed/xXPI4LLrNuM" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
    <br/>


Using pylustrator is very easy and does not require substantial modifications to your code. Just add

.. code-block:: python
    :linenos:

    import pylustrator
    pylustrator.start()

before creating your first figure in your code. When calling ``plt.show()`` the plot will be displayed in a pylustrator
window.

You can test pylustrator with the following example code `example_pylustrator.py <https://raw.githubusercontent.com/rgerum/pylustrator/master/docs/example_pylustrator.py>`_:

.. literalinclude:: example_pylustrator.py
   :language: python
   :emphasize-lines: 6,9
   :linenos:

Saving by pressing ``Ctrl+S`` or confirming to save when closing the window will add some lines of code at the end of your
python script (before your ``plt.show()``) that defines these changes:

.. code-block:: python
    :linenos:

    #% start: automatic generated code from pylustrator
    plt.figure(1).set_size_inches(8.000000/2.54, 8.000000/2.54, forward=True)
    plt.figure(1).axes[0].set_position([0.191879, 0.148168, 0.798133, 0.742010])
    plt.figure(1).axes[0].set_xlabel("data x")
    plt.figure(1).axes[0].set_ylabel("data y")
    plt.figure(1).axes[1].set_position([0.375743, 0.603616, 0.339534, 0.248372])
    plt.figure(1).axes[1].set_xlabel("data x")
    plt.figure(1).axes[1].set_ylabel("data y")
    plt.figure(1).axes[1].set_ylim(-40.0, 90.0)
    #% end: automatic generated code from pylustrator

.. note::
   Because pylustrator can optionally save changes you've made in the GUI to update your source
   code, it cannot be used from a shell. To use pylustrator, call it directly from a
   python file and use the command line to execute.

The good thing is that this generated code is plain matplotlib code, so it will still work when you remove pylustrator
from your code! This is especially useful if you want to distribute your code and do not want to require pylustrator as
a dependency.

Can styling plots be any easier?

Note
----

If you encounter any bugs or unexpected behaviour, you are encouraged to report a bug in our
Github `bugtracker <https://github.com/rgerum/pylustrator/issues>`_.


Citing Pylustrator
------------------

If you use Pylustrator for your publications I would highly appreciate it if you cite the Pylustrator:

* Gerum, R., (2020). **pylustrator: code generation for reproducible figures for publication**. Journal of Open Source Software, 5(51), 1989. `doi:10.21105/joss.01989 <https://doi.org/10.21105/joss.01989>`_


License
-------

Pylustrator is released under the `GPLv3 <https://choosealicense.com/licenses/gpl-3.0/>`_ license. The generated output
code of Pylustrator can be freely used according to the `MIT <https://choosealicense.com/licenses/mit/>`_ license, but as
it relys on Matplotlib also the `Matplotlib License <https://matplotlib.org/users/license.html>`_ has to be taken into
account.

.. toctree::
   :caption: Contents
   :maxdepth: 2

   styling
   composing
   api
API
===

The API of pylustrator is kept quite simple. Most interaction with the pylustrator package is by using the interactive
interface.

.. autofunction:: pylustrator.start

.. autofunction:: pylustrator.load
