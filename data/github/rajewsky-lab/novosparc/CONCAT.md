|PyPI| |Docs| |PePy|

.. |PyPI| image:: https://img.shields.io/pypi/v/novosparc.svg
   :target: https://pypi.org/project/novosparc/
.. |Docs| image:: https://readthedocs.org/projects/novosparc/badge/?version=latest
   :target: https://novosparc.readthedocs.io/
.. |PePy| image:: https://static.pepy.tech/badge/novosparc
   :target: https://pepy.tech/project/novosparc

novoSpaRc - *de novo* Spatial Reconstruction of Single-Cell Gene Expression
===========================================================================

.. image:: https://raw.githubusercontent.com/nukappa/nukappa.github.io/master/images/novosparc.png
   :width: 90px
   :align: left

``novoSpaRc`` predicts locations of single cells in space by solely using 
single-cell RNA sequencing data. An existing reference database of marker genes
is not required, but significantly enhances performance if available.

``novoSpaRc`` accompanies the following publications:

    | *Gene Expression Cartography*
    | M Nitzan*, N Karaiskos*, N Friedman†, N Rajewsky†
    | `Nature (2019) <https://www.nature.com/articles/s41586-019-1773-3>`_

and

    | *novoSpaRc: flexible spatial reconstruction of single-cell gene expression with optimal transport*
    | N Moriel*, E Senel*, N Friedman, N Rajewsky, N Karaiskos†, M Nitzan†
    | `Nature Protocols (2021) <https://www.nature.com/articles/s41596-021-00573-7>`_

Read the `documentation <https://novosparc.readthedocs.io>`_ and the 
`tutorial <https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_drosophila_embryo_tutorial.ipynb>`_ for more information.
.. role:: small
.. role:: smaller
.. role:: noteversion

Version 0.4.4 :small:`11 October 2021`
---------------------------------------
Fixed a bug regarding the distance metric usage in cost calculation.

Version 0.4.3 :small:`20 April 2021`
---------------------------------------
Fixed bugs. Added self consistency analysis and updated tutorials.

Version 0.4.2 :small:`03 April 2021`
---------------------------------------
Improved package structure, fixed minor performace issues, and fixed bugs. Added two new tutorials (corti & osteosarcoma) and 
updated the previous tutorials with validation analyses.

Version 0.4.1 :small:`24 August 2020`
---------------------------------------
Changed the package structure and run flow of the scripts. Added anndata and scanpy support. Updated tutorials and implemented
basic target geometries.

Version 0.3.11 :small:`27 April 2020`
---------------------------------------
Moran's I algorithm for spatially informative genes is implemented and removed pysal dependency.

Version 0.3.10 :small:`07 February 2020`
---------------------------------------
Added Moran's I algorithm to detect spatially informative genes.

Version 0.3.7 :small:`29 October 2019`
---------------------------------------
Updated computation of shortest paths that singificantly reduces
running time.

Version 0.3.5 :small:`13 June 2019`
---------------------------------------
Fixed a bug that was prone to produce infinities during reconstruction.
Improved plotting functions and added new ones for plotting mapped cells.

Version 0.3.4 :small:`27 February 2019`
---------------------------------------
novoSpaRc reconstructs single-cell gene expression without relying on existing
reference markers and makes great use of such information if available.
General usage 
=============
To spatially reconstruct gene expression, ``novoSpaRc`` performs the following
steps:

1. Read the gene expression matrix.
   
   a. *Optional*: select a random set of cells for the reconstruction.
   b. *Optional*: select a small set of genes (e.g. highly variable).
2. Construct the target space.
3. Setup the optimal transport reconstruction.

   a. *Optional*: use existing information of marker genes, if available.
4. Perform the spatial reconstruction.

   a. Assign cells a probability distribution over the target space.
   b. Derive a virtual in situ hybridization (vISH) for all genes over the target space.

5. Write outputs to file for further use, such as the spatial gene expression matrix and the target space coordinates.
6. *Optional*: plot spatial gene expression patterns.
7. *Optional*: identify and plot spatial archetypes.

Demonstration
~~~~~~~~~~~~~
We provide scripts that spatially reconstruct two of the tissues presented
in the paper: the intestinal epithelium [Moor18]_ and the stage 6 Drosophila embryo
[BDTNP]_. 

See also our `tutorial <https://github.com/rajewsky-lab/novosparc/blob/master/reconstruct_drosophila_embryo_tutorial.ipynb>`_ on reconstructing the Drosophila embryo.

The intestinal epithelium
~~~~~~~~~~~~~~~~~~~~~~~~~
The ``reconstruct_intestine_denovo.py`` script reconstructs the crypt-to-villus axis of the mammalian intestinal epithelium, based on data from [Moor18]_. 
The reconstruction is performed *de novo*, without using any marker genes. 
The script outputs plots of (a) a histogram showing the distribution of assignment values over embedded zones for each original villus zone, and (b) average spatial gene expression over the original villus zones and embedded zones of 4 gene groups.

Running time on a standard computer is under a minute.

The *Drosophila* embryo
~~~~~~~~~~~~~~~~~~~~~~~
The ``reconstruct_bdtnp_with_markers.py`` script reconstructs the early
*Drosophila* embryo with only a handful of markers, based on the [BDTNP]_ dataset. 
All cells are used and
a random set of 1-4 markers is selected. The script outputs plots of
gene expression for a list of genes, as well as Pearson correlations of the
reconstructed and original expression values for all genes.
Notice that the results depend on which marker genes are selected. 
In the manuscript we averaged the results over many different choices of marker genes.

Running time on a standard desktop computer is around 6-7 minutes.

Running novoSpaRc on your data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A template file for running ``novoSpaRc`` on custom datasets is 
provided (``reconstruct_tissue.py``). To successfully run ``novoSpaRc`` modify the
template file accordingly.

Constructing different grid shapes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We advise to use ```novoSpaRc`` with diverse target spaces to assess how robust
the spatial reconstructions are. A straightforward way to create a target space
which is more interesting than a square grid, is to have a simple image with the
target space painted in black on it, such as the one below:

.. image:: https://raw.githubusercontent.com/nukappa/nukappa.github.io/master/images/tissue_example.png
   :width: 200px
   :align: center

Then use the function ``create_target_space_from_image`` from the geometry module
to read the image and create a target space out of it. It is advisable to
sample a number of all the read locations and not use them all.

Installation
------------

A working ``Python 3.5`` installation and the following libraries are required: 
``matplotlib``, ``numpy``, ``sklearn``, ``scipy``, ``ot`` and ``networkx``.

The code is partially based on adjustments of the `POT (Python Optimal Transport) <https://github.com/rflamary/POT>`_ library.


``novoSpaRc`` requires a working ``Python 3.4`` installation.


PyPI
~~~~

To install ``novoSpaRc`` try::

    pip install novoSpaRc


Trouble shooting
~~~~~~~~~~~~~~~

If you do not have sudo rights (you get a ``Permission denied`` error)::

    pip install --user novosparc

If installation through ``pip`` fails try installing the ``pot`` library
first::

    pip install cython
    pip install pot

and then ``novoSpaRc``::

    pip install novosparc 
.. novosparc documentation master file, created by
   sphinx-quickstart on Fri Feb  8 11:23:52 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: ../README.rst
   :end-line: 24

.. include:: release_notes.rst

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   tutorials
   references
References
----------

.. [BDTNP] BDTNP,
   *Berkeley Drosophila Transcription Network Project*,
   `<bdtnp.lbl.gov>`__.

.. [Halpern17] Halpern *et al.* (2017),
   *Single-cell spatial reconstruction reveals global division of labour in the mammalian liver*,
   `Nature <https://doi.org/10.1038/nature21065>`__.

.. [Karaiskos17] Karaiskos *et al.* (2017),
   *The Drosophila embryo at single-cell transcriptome resolution*,
   `Science <https://doi.org/10.1126/science.aan3235>`__.

.. [Moor18] Moor *et al.* (2018),
   *Spatial Reconstruction of Single Enterocytes Uncovers Broad Zonation along the Intestinal Villus Axis*,
   `Cell <https://doi.org/10.1016/j.cell.2018.08.063>`__.

.. [Nitzan18] Nitzan *et al.* (2018),
   *Charting tissues from single-cell transcriptomes*,
   `bioRxiv <https://www.biorxiv.org/content/10.1101/456350v1>`__.
