mpi4jax
=======

|JOSS paper| |PyPI Version| |Conda Version| |Tests| |codecov| |Documentation Status|

``mpi4jax`` enables zero-copy, multi-host communication of `JAX <https://jax.readthedocs.io/>`_ arrays, even from jitted code and from GPU memory.


But why?
--------

The JAX framework `has great performance for scientific computing workloads <https://github.com/dionhaefner/pyhpc-benchmarks>`_, but its `multi-host capabilities <https://jax.readthedocs.io/en/latest/jax.html#jax.pmap>`_ are still limited.

With ``mpi4jax``, you can scale your JAX-based simulations to *entire CPU and GPU clusters* (without ever leaving ``jax.jit``).

In the spirit of differentiable programming, ``mpi4jax`` also supports differentiating through some MPI operations.


Installation
------------

``mpi4jax`` is available through ``pip`` and ``conda``:

.. code:: bash

   $ pip install mpi4jax                     # Pip
   $ conda install -c conda-forge mpi4jax    # conda

If you use pip and don't have JAX installed already, you will also need to do:

.. code:: bash

   $ pip install jaxlib

(or an equivalent GPU-enabled version, `see the JAX installation instructions <https://github.com/google/jax#installation>`_)

In case your MPI installation is not detected correctly, `it can help to install mpi4py separately <https://mpi4py.readthedocs.io/en/stable/install.html>`_. When using a pre-installed ``mpi4py``, you *must* use ``--no-build-isolation`` when installing ``mpi4jax``:

.. code:: bash

   # if mpi4py is already installed
   $ pip install cython
   $ pip install mpi4jax --no-build-isolation

`Our documentation includes some more advanced installation examples. <https://mpi4jax.readthedocs.io/en/latest/installation.html>`_


Example usage
-------------

.. code:: python

   from mpi4py import MPI
   import jax
   import jax.numpy as jnp
   import mpi4jax

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()

   @jax.jit
   def foo(arr):
      arr = arr + rank
      arr_sum, _ = mpi4jax.allreduce(arr, op=MPI.SUM, comm=comm)
      return arr_sum

   a = jnp.zeros((3, 3))
   result = foo(a)

   if rank == 0:
      print(result)

Running this script on 4 processes gives:

.. code:: bash

   $ mpirun -n 4 python example.py
   [[6. 6. 6.]
    [6. 6. 6.]
    [6. 6. 6.]]

``allreduce`` is just one example of the MPI primitives you can use. `See all supported operations here <https://mpi4jax.readthedocs.org/en/latest/api.html>`_.

Community guidelines
--------------------

If you have a question or feature request, or want to report a bug, feel free to `open an issue <https://github.com/mpi4jax/mpi4jax/issues>`_.

We welcome contributions of any kind `through pull requests <https://github.com/mpi4jax/mpi4jax/pulls>`_. For information on running our tests, debugging, and contribution guidelines please `refer to the corresponding documentation page <https://mpi4jax.readthedocs.org/en/latest/developers.html>`_.

How to cite
-----------

If you use ``mpi4jax`` in your work, please consider citing the following article:

::

  @article{mpi4jax,
    doi = {10.21105/joss.03419},
    url = {https://doi.org/10.21105/joss.03419},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {65},
    pages = {3419},
    author = {Dion Häfner and Filippo Vicentini},
    title = {mpi4jax: Zero-copy MPI communication of JAX arrays},
    journal = {Journal of Open Source Software}
  }

.. |Tests| image:: https://github.com/mpi4jax/mpi4jax/workflows/Tests/badge.svg
   :target: https://github.com/mpi4jax/mpi4jax/actions?query=branch%3Amaster
.. |codecov| image:: https://codecov.io/gh/mpi4jax/mpi4jax/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/mpi4jax/mpi4jax
.. |PyPI Version| image:: https://img.shields.io/pypi/v/mpi4jax
   :target: https://pypi.org/project/mpi4jax/
.. |Conda Version| image:: https://img.shields.io/conda/vn/conda-forge/mpi4jax.svg
   :target: https://anaconda.org/conda-forge/mpi4jax
.. |Documentation Status| image:: https://readthedocs.org/projects/mpi4jax/badge/?version=latest
   :target: https://mpi4jax.readthedocs.io/en/latest/?badge=latest
.. |JOSS paper| image:: https://joss.theoj.org/papers/10.21105/joss.03419/status.svg
   :target: https://doi.org/10.21105/joss.03419
Demo application: Shallow-water model
=====================================

To show you what ``mpi4jax`` is capable of, we include a full implementation of a physical `nonlinear shallow-water model <https://github.com/dionhaefner/shallow-water>`_.

A shallow-water model simulates the evolution of the sea surface if temperature and salinity of the water do not vary with depth. Our nonlinear implementation is even capable of modelling turbulence. A possible solution looks like this:

.. raw:: html

    <video width="80%" style="margin: 1em auto; display: block;" controls>
        <source src="_static/shallow-water.mp4" type="video/mp4">
        Your browser does not support the video tag.
    </video>


The demo script is too long to include here, but you can
:download:`download it <../examples/shallow_water.py>` or :doc:`see the source here <shallow-water-source>`.


Running the demo
----------------

Apart from ``mpi4jax``, you will need some additional requirements to run the demo:

.. code:: bash

    $ pip install matploblib tqdm

Then, you can run it like this:

.. code:: bash

    $ mpirun -n 4 python shallow_water.py
    WARNING:absl:No GPU/TPU found, falling back to CPU. (Set TF_CPP_MIN_LOG_LEVEL=0 and rerun for more info.)
    WARNING:absl:No GPU/TPU found, falling back to CPU. (Set TF_CPP_MIN_LOG_LEVEL=0 and rerun for more info.)
    WARNING:absl:No GPU/TPU found, falling back to CPU. (Set TF_CPP_MIN_LOG_LEVEL=0 and rerun for more info.)
    WARNING:absl:No GPU/TPU found, falling back to CPU. (Set TF_CPP_MIN_LOG_LEVEL=0 and rerun for more info.)
    100%|█████████▉| 9.98/10.00 [00:28<00:00,  2.90s/model day]
    Solution took 25.79s

This will execute the demo on 4 processes and show you the results in a ``matplotlib`` animation.


Benchmarks
----------

Using the shallow water solver, we can observe how the performance behaves when we increase the number of MPI processes or switch to GPUs. Here we show some benchmark results on a machine with 2x Intel Xeon E5-2650 v4 CPUs and 2x NVIDIA Tesla P100 GPUs.

.. note::

    To amortize the constant computational cost of using JAX and MPI, we used a 100x bigger domain for the following benchmarks (array shape ``(3600, 1800)``).

.. code:: bash

    # CPU
    $ JAX_PLATFORM_NAME=cpu mpirun -n 1 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [01:55<00:09, 1248.13s/model day]
    Solution took 111.95s

    $ JAX_PLATFORM_NAME=cpu mpirun -n 2 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [01:33<00:07, 1010.01s/model day]
    Solution took 89.67s

    $ JAX_PLATFORM_NAME=cpu mpirun -n 4 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [00:41<00:03, 451.75s/model day]
    Solution took 38.57s

    $ JAX_PLATFORM_NAME=cpu mpirun -n 6 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [00:31<00:02, 345.56s/model day]
    Solution took 28.70s

    $ JAX_PLATFORM_NAME=cpu mpirun -n 8 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [00:23<00:01, 260.17s/model day]
    Solution took 20.62s

    $ JAX_PLATFORM_NAME=cpu mpirun -n 16 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [00:19<00:01, 208.55s/model day]
    Solution took 15.73s

    # GPU
    $ JAX_PLATFORM_NAME=gpu mpirun -n 1 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [00:09<00:00, 103.18s/model day]
    Solution took 6.28s

    $ JAX_PLATFORM_NAME=gpu MPI4JAX_USE_CUDA_MPI=0 mpirun -n 2 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [00:07<00:00, 76.42s/model day]
    Solution took 3.87s

    $ JAX_PLATFORM_NAME=gpu MPI4JAX_USE_CUDA_MPI=1 mpirun -n 2 -- python examples/shallow_water.py --benchmark
    92%|█████████▏| 0.09/0.10 [00:07<00:00, 76.28s/model day]
    Solution took 3.89s


Using pure NumPy, the same model takes over 12 minutes (770s) to execute on 1 process, which goes to show how efficient JAX is for this kind of workload.
Developer guide
===============

Development install
-------------------

To install mpi4jax in editable mode along with all optional dependencies for testing, just run

.. code:: bash

   $ pip install -e .[dev]

from the repository root.

Running tests
-------------

We use ``pytest`` for testing. After installing the development dependencies, you can run our testing suite with the following commands:

.. code:: bash

    $ pytest .
    $ mpirun -np 2 pytest .

.. warning::

    Just executing ``pytest`` will run the tests on only 1 process, which means that a large part of mpi4jax cannot be tested (because it relies on communication between different processes). Therefore, you should always make sure that the tests also pass on multiple processes (via ``mpirun``).

Contributing
------------

We welcome code contributions or changes to the documentation via `pull requests <https://github.com/mpi4jax/mpi4jax/pulls>`_ (PRs).

We use pre-commit hooks to enforce a common code format. To install
them, just run:

.. code:: bash

   $ pre-commit install

in the repository root. Then, all changes will be validated automatically before you commit.

.. note::

    If you introduce new code, please make sure that it is covered by tests.
    To catch problems early, we recommend that you run the test suite locally before creating a PR.

Debugging
---------

You can set the environment variable ``MPI4JAX_DEBUG`` to ``1`` to
enable debug logging every time an MPI primitive is called from within a
jitted function. You will then see messages like these:

.. code:: bash

   $ MPI4JAX_DEBUG=1 mpirun -n 2 python send_recv.py
   r0 | MPI_Send -> 1 with tag 0 and token 7fd7abc5f5c0
   r1 | MPI_Recv <- 0 with tag -1 and token 7f9af7419ac0

This can be useful to debug deadlocks or MPI errors.
Usage examples
==============

Basic example: Global sum
-------------------------

The following computes the sum of an array over several processes (similar to :func:`jax.lax.psum`), using :func:`~mpi4jax.allreduce`:

.. code:: python

   from mpi4py import MPI
   import jax
   import jax.numpy as jnp
   import mpi4jax

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()

   @jax.jit
   def foo(arr):
      arr = arr + rank
      arr_sum, _ = mpi4jax.allreduce(arr, op=MPI.SUM, comm=comm)
      return arr_sum

   a = jnp.zeros((3, 3))
   result = foo(a)

   if rank == 0:
      print(result)

Most MPI libraries supply a wrapper executable ``mpirun`` to execute a script on several processes:

.. code:: bash

   $ mpirun -n 4 python mpi4jax-example.py
   [[6. 6. 6.]
    [6. 6. 6.]
    [6. 6. 6.]]

The result is an array full of the value 6, because each process adds its rank to the result (4 processes with ranks 0, 1, 2, 3).

Basic example: sending and receiving
------------------------------------

``mpi4jax`` can of course also send and receive data without performing an operation on it. For this, you can use :func:`~mpi4jax.send` and :func:`~mpi4jax.recv`:

.. _example_2:

.. code:: python

    from mpi4py import MPI
    import jax
    import jax.numpy as jnp
    import mpi4jax

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    assert size == 2  # make sure we are on 2 processes

    @jax.jit
    def foo(arr):
        arr = arr + rank
        # note: this could also use mpi4jax.sendrecv
        if rank == 0:
            # send, then receive
            token = mpi4jax.send(arr, dest=1, comm=comm)
            other_arr, token = mpi4jax.recv(arr, source=1, comm=comm, token=token)
        else:
            # receive, then send
            other_arr, token = mpi4jax.recv(arr, source=0, comm=comm)
            token = mpi4jax.send(arr, dest=0, comm=comm, token=token)

        return other_arr

    a = jnp.zeros((3, 3))
    result = foo(a)

    print(f'r{rank} | {result}')

Executing this shows that each process has received the data from the other process:

.. code:: bash

    $ mpirun -n 2 python mpi4jax-example-2.py
    r1 | [[0. 0. 0.]
     [0. 0. 0.]
     [0. 0. 0.]]
    r0 | [[1. 1. 1.]
     [1. 1. 1.]
     [1. 1. 1.]]

For operations like this, the correct order of the :func:`~mpi4jax.send` / :func:`~mpi4jax.recv` calls is critical to prevent the program from deadlocking (e.g. when both processes wait for data at the same time).

In ``mpi4jax``, we enforce order of execution through *tokens*. In :ref:`the example code <example_2>`, you can see this behavior e.g. in the following lines:

.. code:: python

    token = mpi4jax.send(arr, dest=1, comm=comm)
    other_arr, token = mpi4jax.recv(arr, source=1, comm=comm, token=token)

The first call to :func:`~mpi4jax.send` returns a token, which we then pass to :func:`~mpi4jax.recv`. :func:`~mpi4jax.recv` *also* returns a new token that we could pass to subsequent communication primitives.

Because of the nature of JAX, **using tokens to enforce order is not optional.** If you do not use correct token management, you will experience deadlocks and crashes.

.. seealso::

    For more information on tokens, see :ref:`tokens`.
Installation
============

Basic installation
------------------

Start by `installing a suitable version of JAX and jaxlib <https://github.com/google/jax#installation>`_. If you don't plan on using ``mpi4jax`` on GPU, the following will do:

.. code:: bash

   $ pip install jax jaxlib

.. note::

   Much of the functionality we need has recently been added to JAX, which itself is changing frequently. Therefore, ``mpi4jax`` has somewhat strict requirements on the supported versions of JAX and jaxlib. Be prepared to upgrade!

We recommend that you use ``pip`` to install mpi4jax (but a distribution is also available via ``conda`` which will work if MPI, mpi4py and mpi4jax are all installed through conda ):

.. code:: bash

   $ pip install mpi4jax
   $ conda install -c conda-forge mpi4jax

Installing via ``pip`` requires a working installation of MPI to succeed. If you don't already have MPI and want to get started as quickly as possible, try ``conda``, which bundles the MPI library (but remember `not to mix pip and conda <https://www.anaconda.com/blog/using-pip-in-a-conda-environment>`_).

.. warning::

   We advise against using the conda installation in HPC environments because it is not possible to change the MPI library ``mpi4py`` is linked against.

And that is it! If you are familiar with MPI, you should in principle be able to get started right away. However, we recommend that you have a look at :doc:`sharp-bits`, to make sure that you are aware of some of the pitfalls of ``mpi4jax``.

Selecting the MPI distribution
------------------------------

``mpi4jax`` will use the MPI distribution with which ``mpi4py`` was built.
If ``mpi4py`` is not installed, it will be installed automatically before
installing ``mpi4jax``.

.. warning::

   If ``mpi4py`` is already installed, you *must* use ``--no-build-isolation`` when installing ``mpi4jax``:

   .. code:: bash

      # if mpi4py is already installed
      $ pip install cython
      $ pip install mpi4jax --no-build-isolation

To check which MPI library both libraries link to, run the following command in your
prompt.

.. code:: bash

	$ python -c "import mpi4py; print(mpi4py.get_config())"

If you wish to use a specific MPI library (only possible when using ``pip``), it is
usually sufficient to specify the ``MPICC`` environment variable *before* installing
``mpi4py``.

.. seealso::

   In doubt, please refer to `the mpi4py documentation <https://mpi4py.readthedocs.io/en/stable/install.html>`_.


Installation with GPU support
-----------------------------

.. note::

   To use JAX on the GPU, make sure that your ``jaxlib`` is `built with CUDA support <https://github.com/google/jax#installation>`_.

``mpi4jax`` also supports JAX arrays stored in GPU memory.

To build ``mpi4jax``'s GPU extensions, we need to be able to locate the CUDA headers on your system. If they are not detected automatically, you can set the environment variable :envvar:`CUDA_ROOT` when installing ``mpi4jax``::

   $ CUDA_ROOT=/usr/local/cuda pip install mpi4jax

This is sufficient for most situations. However, ``mpi4jax`` will copy all data from GPU to CPU and back before and after invoking MPI.

If this is a bottleneck in your application, you can build MPI with CUDA support and *communicate directly from GPU memory*. This requires that you re-build the entire stack:

- Your MPI library, e.g. `OpenMPI <https://www.open-mpi.org/faq/?category=buildcuda>`_, with CUDA support.
- ``mpi4py``, linked to your CUDA-enabled MPI installation.
- ``mpi4jax``, using the correct ``mpi4py`` installation.

.. seealso::

   Read :ref:`here <gpu-usage>` on how to use zero-copy GPU communication after installation.
.. include:: ../README.rst

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   usage
   shallow-water
   sharp-bits
   api
   developers

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
🔪 The Sharp Bits 🔪
====================

Read ahead for some pitfalls, counter-intuitive behavior, and sharp edges that we had to introduce in order to make this work.

.. _tokens:

Token management
----------------

The compiler behind JAX, XLA, is not aware of the fact that MPI function calls such as :func:`~mpi4jax.send` or :func:`~mpi4jax.recv` must be performed in a specific order (in jargon, that they have *side-effects*). It is therefore free to reorder those calls. Reordering of MPI calls usually leads to deadlocks, e.g. when both processes end up receiving before sending (instead of send-receive, receive-send).

*Tokens* are JAX's way to ensure that XLA does not re-order statements with side effects by injecting a fake data dependency between them.

This means that you *have* to use proper token management to prevent reordering from occurring. Every communication primitive in ``mpi4jax`` returns a token as the last return object, which you have to pass to subsequent primitives within the same JIT block, like this:

.. code:: python

    # DO NOT DO THIS
    mpi4jax.send(arr, comm=comm)
    new_arr, _ = mpi4jax.recv(arr, comm=comm)

    # INSTEAD, DO THIS
    token = mpi4jax.send(arr, comm=comm)
    new_arr, token = mpi4jax.recv(arr, comm=comm, token=token)

Those functions will then be executed in the same order as the sequence of tokens, from first to last.


No in-place operations in JAX
-----------------------------

JAX arrays are immutable, which means that functions cannot modify their input arguments. Therefore, unlike in ``mpi4py``, operations like :func:`mpi4jax.recv` use their first argument only to determine the correct shape and dtype of the output, but do not populate it with data.

This means that you *cannot* do:

.. code:: python

    # DO NOT DO THIS
    recv_arr = jnp.zeros((10, 10))
    mpi4jax.recv(recv_arr, comm=comm)
    # recv_arr still only contains 0

Instead, you need to use the returned array from :func:`mpi4jax.recv`:

.. code:: python

    # INSTEAD, DO THIS
    recv_arr = jnp.zeros((10, 10))
    recv_arr, _ = mpi4jax.recv(recv_arr, comm=comm)

.. _gpu-usage:

Using CUDA MPI
--------------

``mpi4jax`` is able to communicate data directly from and to GPU memory. :doc:`This requires that MPI, JAX, and mpi4jax are built with CUDA support. <installation>`

Currently, we cannot detect whether MPI was built with CUDA support.
Therefore, by default, ``mpi4jax`` will not read directly from GPU
memory, but instead copy to the CPU and back.

If you are certain that the underlying MPI library was built with CUDA
support, you can set the following environment variable:

.. code:: bash

   $ export MPI4JAX_USE_CUDA_MPI=1

Data will then be copied directly from GPU to GPU. If your MPI library
does not have CUDA support, you will receive a segmentation fault when
trying to access GPU memory.


Using ``mpi4jax`` *and* ``mpi4py``
----------------------------------

.. warning::

    Do not use ``mpi4jax`` and ``mpi4py`` with the same communicator!

Consider the following example, where one process sends some Python data via ``mpi4py`` and JAX data via ``mpi4jax``, and the other process receives it:

.. code:: python

    # DO NOT DO THIS
    import numpy as np
    import jax.numpy as jnp

    from mpi4py import MPI
    import mpi4jax

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    arr_np = np.random.rand(10, 10)
    arr_jax = jnp.zeros((10, 10))

    if rank == 0:
        mpi4jax.send(arr_jax, comm=comm)
        comm.send(arr_np)
    else:
        arr_jax = mpi4jax.recv(arr_jax, comm=comm)
        arr = comm.recv(arr_np)

Because everything is lazily executed in JAX, we cannot rely on a particular execution order. Specifically, we don't know whether the function ``mpi4jax.send`` wille be executed before or after the ``comm.send`` call. In the worst case, this creates a deadlock.

The simplest solution is therefore to stick to *either* ``mpi4py`` *or* ``mpi4jax``. But if you have to use both, make sure that they use different communicators:


.. code:: python

    # INSTEAD, DO THIS
    import numpy as np
    import jax.numpy as jnp

    from mpi4py import MPI
    import mpi4jax

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # create a new communicator for mpi4jax
    comm_jax = comm.Clone()

    arr_np = np.random.rand(10, 10)
    arr_jax = jnp.zeros((10, 10))

    if rank == 0:
        mpi4jax.send(arr_jax, comm=comm_jax)
        comm.send(arr_np)
    else:
        arr_jax = mpi4jax.recv(arr_jax, comm=comm_jax)
        arr = comm.recv(arr_np)

    comm_jax.Free()
:orphan:

:file:`examples/shallow_water.py`
=================================

.. literalinclude:: ../examples/shallow_water.py
   :language: python
   :linenos:

:download:`Download source (shallow_water.py) <../examples/shallow_water.py>`
Utilities
=========

has_cuda_support
----------------

.. autofunction:: mpi4jax.has_cuda_support


Communication primitives
========================

allgather
---------

.. autofunction:: mpi4jax.allgather

allreduce
---------

.. autofunction:: mpi4jax.allreduce

alltoall
--------

.. autofunction:: mpi4jax.alltoall

barrier
-------

.. autofunction:: mpi4jax.barrier

bcast
-----

.. autofunction:: mpi4jax.bcast

gather
------

.. autofunction:: mpi4jax.gather

recv
----

.. autofunction:: mpi4jax.recv

reduce
------

.. autofunction:: mpi4jax.reduce

scan
----

.. autofunction:: mpi4jax.scan

scatter
-------

.. autofunction:: mpi4jax.scatter

send
----

.. autofunction:: mpi4jax.send

sendrecv
--------

.. autofunction:: mpi4jax.sendrecv
