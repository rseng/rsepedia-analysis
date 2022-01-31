# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
### Added
- Empty Python project directory structure



# PSRDada Python

[![DOI](https://zenodo.org/badge/132437444.svg)](https://zenodo.org/badge/latestdoi/132437444)

Python3 bindings to the ringbuffer implementation in [PSRDada](http://psrdada.sourceforge.net/)

It allows you to connect to a ringbuffer, read header info, and read/write ringbuffer pages as numpy arrays.
I aim to be interoperable, and bug-for-bug compattible, with PSRDADA.
It also supports reading multiple datasets from a ringbuffer for connecting to dada\_dbevent.

This is a proof-of-concept implementation, only functions from PSRDADA immediately necessary for my project will be implemented.
The code is inspired by the example code in the [cython userguide](http://cython.readthedocs.io/en/latest/src/userguide/buffer.html)

Error messages, together with a start-up and shutdown message, are written to the syslog.
Depending on your OS, you can read those using:

```bash
 $ grep psrdada-python /var/log/syslog
 $ journalctl --since "10 minutes ago" | grep psrdada-python
```

# Documentation

[The documentation is available here](https://trasal.github.io/psrdada-python/)
See the two scripts in the examples directory, and the tests.

```python
#!/usr/bin/env python3
import numpy as np
from psrdada import Reader

# Create a reader instace and connect to a running ringbuffer with key 'dada'
reader = Reader(0xdada)

# loop over the pages
for page in reader:
    # read the page as numpy array
    data = np.asarray(page)
    print(np.sum(data))

reader.disconnect()
```

# Installation
clone the repository  
    `git clone git@github.com:TRASAL/psrdada-python.git`
change into the top-level directory  
    `cd psrdada-python`  
install the dependencies (in a virtual env)
    `python -m venv env && . env/bin/activate && pip install -f requirements.txt`
build the package
    `make && make test && make install`

# Dependencies

PSRDada, see their [website](https://sourceforge.net):

 * PSRDada dada\_db exectuable in the PATH for testing;

 * PSRDada header files needed for compilation, set CPATH or CFLAGS.

 * PSRDada library needed during runtime, set LD\_LIBRARY\_PATH

# License
Copyright (c) 2018, Jisk Attema
Apache Software License 2.0.

# Contributing

All contributions are welcome!
Please use the github issue tracker to get in touch.

Contributing authors so far:
* Jisk Attema
* Leon Oostrum
* Liam Connor
* Wael Farah

Contributing
============

To contribute to PSRDada Python please contact Jisk Attema <j.attema@esciencecenter.nl>


Index
=====

`genindex`
.. psrdada-python documentation master file, created by
   sphinx-quickstart on Wed Jun  5 16:04:58 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Api documentation
=================

.. automodule:: psrdada


Ringbuffer
----------

.. py:class:: Ringbuffer

   .. py:attribute:: isConnected

      Boolean; indicates if the Ringbuffer instance is connected to a
      running dada_db.

   .. py:attribute:: isHoldingPage

      Boolean; indicates if the Ringbuffer instance is currently holding a page;
      ie. has exclusive write access to it for Writers, or is preventing it from 
      being rewritten for Readers.

   .. py:attribute:: isEndOfData

      Boolean; indicates if the EOD flag has been set on the Ringbuffer. Used
      to implement iterators for the Reader and Writer classes.

   .. py:attribute:: header

      dict; contains a copy of the last read header. For Readers, the original
      un-parsed header is available under the __RAW_HEADER__ key.

Reader
------

.. autoclass:: Reader
   :show-inheritance:
   :members:

Writer
------

.. autoclass:: Writer
   :show-inheritance:
   :members:

Viewer
------

.. autoclass:: Viewer
   :show-inheritance:
   :members:
.. psrdada-python documentation master file, created by
   sphinx-quickstart on Wed Jun  5 16:04:58 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to psrdada-python's documentation
=========================================

.. automodule:: psrdada
   :noindex:

.. toctree::

   basic
   apidoc
   genindex
Basic example
-------------

Using iterators
~~~~~~~~~~~~~~~

Connect to and read from a PSRDada ringbuffer in a pythonic way, using
the Iterator interface:

.. literalinclude:: ../examples/sum_with_iterator.py

Explicit ringbuffer control 
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For usecases that require a more finegrained flow control, there is a
more C-like API:

.. literalinclude:: ../examples/sum_no_iterator.py

More examples
~~~~~~~~~~~~~

For more examples, see the *examples* and *tests* directories in the repo.
Examples cover writing to the ringbuffer, reading and writing header data, and dealing with EndOfData.
EndOfData is used by *dada_dbevent* to send independent datasets through a ringbuffer, instead of a single continuous stream.
