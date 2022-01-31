Coding style
============

* We follow `PEP 8 <http://legacy.python.org/dev/peps/pep-0008/>`_.

* Before you check in code, please run ``pyflakes`` and ``pep8``
  (or ``flake8``) on your source files, and make sure you don't get any
  errors.

* Please don't use backticks, they're deprecated. To format strings, use
  the ``%`` operator or ``str.format``.

* Run tests using ``nosetests``.

* Please don't use ``assert`` in tests. Use the ``numpy.testing`` and
  ``nose.tools`` assertions.

* Put docstrings on classes, not initializers (``__init__``). Python tools
  expect docstrings on classes.
|Travis|_ |Quality-score|_ |Coverage|_

.. |Travis| image:: https://api.travis-ci.org/NLeSC/PattyAnalytics.png?branch=master
.. _Travis: https://travis-ci.org/NLeSC/PattyAnalytics

.. |Quality-score| image:: https://scrutinizer-ci.com/g/NLeSC/PattyAnalytics/badges/quality-score.png?b=master
.. _Quality-score: https://scrutinizer-ci.com/g/NLeSC/PattyAnalytics/

.. |Coverage| image:: https://scrutinizer-ci.com/g/NLeSC/PattyAnalytics/badges/coverage.png?b=master
.. _Coverage: https://scrutinizer-ci.com/g/NLeSC/PattyAnalytics/

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.45925.svg
   :target: http://dx.doi.org/10.5281/zenodo.45925

Patty Analytics
===============

Reusable point cloud analytics software. Includes segmentation, registration,
file format conversion. This makes use of the python bindings of the
Point Cloud Library (PCL; <https://github.com/NLeSC/python-pcl>).

Copyright 2014-2015 Netherlands eScience Center. Covered by the Apache 2.0
License, see the file ``LICENSE.txt``.

Installing
----------

First install the required dependencies:

* Python 2.7
* NumPy
* SciPy
* virtualenv
* LibLAS
* PCL 1.7

Now set up an environment::

    $ virtualenv /some/where/env --system-site-packages
    $ . /some/where/env/activate


Install the python packages listed in ``requirements.txt`` using pip::

    $ pip install -r requirements.txt
    $ pip install -U nose  # make sure nose runs in the virtualenv
    $ python setup.py install

To exit the python virtualenv run::

    $ deactivate

Running
-------
The main functionality of PatTy Analytics is contained in the **registration**
script. This script takes an unaligned dense point cloud and attempts to
align in on to the existing drivemap. The script also attempts to find the
optimal scaling, rotation and orientation of the dense point cloud, as part of
the alignment process. This script can be run as follows::

    $ python scripts/registration.py SOURCE.las DRIVEMAP.las FOOTPRINT.csv OUTPUT.las

where:

  - *SOURCE.las* -- is the dense point cloud to be registered
  - *DRIVEMAP.las* -- is the drive map where the point cloud is registered to.
  - *FOOTPRINT.csv* -- is the footprint of the point cloud on the drive map.
  - *OUTPUT.las* -- is the resulting registered point cloud.

additionally, an *upfile.json* containing the up vector (estimated from the
camera position) can be provided.

    $ python scripts/registration.py SOURCE.las DRIVEMAP.las FOOTPRINT.csv OUTPUT.las -u upfile.json

## Examples

The following image is a screenshot of a dense point cloud to be registered
on the drive map -- this would correspond to *SOURCE.las*.

![Site 558 dense point cloud](./img/site558_dense.png?raw=true "Dense point cloud")

This screenshot shows the drive map where we want to register to -- this corresponds
to *DRIVEMAP.las*

![Site 558 drive map point cloud](./img/site558_drivemap.png?raw=true "Drive map")

Finally, this screenshot shows the dense point cloud registered on the drive map.
The dense point cloud has been rotated, scaled and translated to find it's best
fit on the drive map -- this corresponds to *OUTPUT.las*.


![Site 558 registered point cloud](./img/site558_registered.png?raw=true "Registered point cloud")

Testing
-------

To run unit tests, issue::

    $ nosetests

Documentation
-------------

Documentation can be found here_

.. _here: http://nlesc.github.io/PattyAnalytics/
.. Patty Analytics documentation master file, created by
   sphinx-quickstart on Wed Mar 18 10:04:17 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Patty Analytics
===============

Patty Analytics aims to register pointclouds that were generated from photos
or video to an absolute position, scale and orientation.

Pointclouds generated from photos are generally messy; they have holes and
floating unidentified objects. In our scripts we assume to have the following
information: a map (drivemap) which has an extremely low resolution but has
good absolute coordinates; a footprint polygon denoting more or less the
latitude and longitude and area of the object (x and y coordinates). Finally,
we have the high-resolution pointcloud of the object. By the nature of
creating this pointcloud, it is densest at the object, since the photos
usually center on this object. In some cases, there are also camera positions
available, relative to the object.


Registration
============

The `registration.py` script is the main registration algorithm. This script can be run as follows:

```
Usage:
  registration.py [-h] [-d <sample>] [-u <upfile>] <source> <drivemap> <footprint> <output>

Positional arguments:
  source       Source LAS file
  drivemap     Target LAS file to map source to
  footprint    Footprint for the source LAS file
  output       file to write output LAS to
```

Example run:
```
python scripts/registration.py data/SOURCE/SITE_X.las data/DRIVEMAP/X.las data/FOOTPRINT/X.footprint.csv data/OUTPUT/SITE_X.las -d 100000
```

Initialization
--------------

The registration script applies a dbscan clustering algorithm before anything
else to filter out noise and end up with the densest parts of the pointcloud.
This is the slowest operation. The algorithm includes at least 70% of the
pointcloud. To speed this operation up and make its runtime predictable, the
pointcloud can be subsampled to a given number of points before doing the
clustering algorithm. A uniform probability distribution is used for this.

Orientation
-----------

The registration script does not yet give a definite orientation. It assumes
that objects stand on the ground, and that the footprint represents some part
of the ground. In turn we assume that both the pointcloud and drivemap are
more or less flat at this location. The first step is then to detect the
boundary of the object, resulting in a line-like subset of the pointcloud,
usually representing its demarcation at the ground and holes that have formed
during pointcloud generation. Then, a small margin around the footprint is cut
out of the drivemap. The orientation is then determined by aligning the two
principal axes of the pointclouds boundary to the drivemap cutout. The
principal axes are found by taking the first two components of a PCA of the
drivemap and the cutout. With these operations, it is still not possible to
distinguish up from down and front from back.

The next step is to, when available, use the camara positions to determine an
approximate UP orientation. This orientation is used to determine up from
down. Unfortunately this still does not determine front from back.

Scale
-----

In our project, red-white meter sticks of 80 cm were placed next to the
objects, so this is a natural way to determine the object's scale. These
sticks have two red and two white parts, each 20 cm long. The red parts are
rare in many natural environments, so we first take all red points (in HSV
color space: H > 0.9 and S > 0.5). Then, the remaining points are clustered
with the dbscan algorithm. The length of each cluster (each red stick part) is
determined by applying principal component analysis to find the longest axis
and measuring this axis. If three or more clusters have about the same length
and high number of points, we assign a high confidence level to the stick
scale. The fewer points we have and the fewer conforming clusters, the lower
our confidence that this scale is the correct one.

If these meter sticks are not available or not well represented, we determine
the scale by comparing the size of the oriented boundary of the object to the
size of the footprint, and take that as the scale.

Position
--------

Once scaled and oriented, the position of the pointcloud is determined by
shifting the pointcloud boundary to the location of the footprint. Its height
is determined by the height of the cutout of the drivemap around the footprint.

Future work
-----------

Once the pointcloud is placed at its approximate position, a few measures
could be taken to improve its position. One way is to use ICP or GICP from the
PCL library. Since the density of the drivemap and the pointcloud differ
greatly, the poincloud probably has to be subsampled first. A complicating
factor is that the drivemap may have a grid-like topology because of its low
density whereas the pointcloud does not, which is likely to throw an ICP
algorithm off.

Additionally, a measure of fit could be used (nearest point cumulative
distance comes to mind) to try a few common operations. The following
operations come to mind: rotating the front to the back; scaling between 50%
and 150%; and slightly moving the object to each side.

The Python-PCL library currently does not support point normals. If it did, it
would give additional tools to determine rotation.

Pointcloud utilities
====================

In addition to registration, this package has a few utilities for dealing with
pointclouds. In `transform.py`, pointclouds can be rotated, scaled and
translated. With `convert.py`, a pointcloud can be converted from one format
(of PCD, PLY or LAS) to another. The Spatial Reference System (SRS) of a LAS
file can be set with `las_set_srs.py`. Point normals are not preserved in
these operations. In the `tests/helpers.py` file, synthetic pointclouds can be
generated for means of testing.

Contents:

.. toctree::
   :maxdepth: 2

   api



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. _api:

API reference
=============

Registration
------------

.. automodule:: patty.registration
    :members:
    :imported-members:

Segmentation
------------

.. automodule:: patty.segmentation
    :members:
    :imported-members:

SRS
---

.. automodule:: patty.srs
    :members:

Utils
-----

.. automodule:: patty.utils
    :members:
