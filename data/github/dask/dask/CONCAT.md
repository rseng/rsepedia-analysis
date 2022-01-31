Dask is a community maintained project. We welcome contributions in the form of bug reports, documentation, code, design proposals, and more. 

For general information on how to contribute see https://docs.dask.org/en/latest/develop.html.
See [developer documentation](https://docs.dask.org/en/latest/develop.html)
for tips on how to get started.
- [ ] Closes #xxxx
- [ ] Tests added / passed
- [ ] Passes `pre-commit run --all-files`
## HDFS testing on Travis CI

Dask & HDFS testing relies on a docker container. The tests are setup to run on
Travis CI, but only under the following conditions:

- Merges to main
- PRs where the commit message contains the string `"test-hdfs"`

If you make a PR changing HDFS functionality it'd be good to have the HDFS
tests run, please add `"test-hdfs"` to your commit message.

## Setting up HDFS testing locally using Docker

Assumes docker is already installed and the docker-daemon is running.

From the root directory in the repo:

- First get the docker container:

```bash
# Either pull it from docker hub
docker pull daskdev/dask-hdfs-testing

# Or build it locally
docker build -t daskdev/dask-hdfs-testing continuous_integration/hdfs/
```

- Start the container and wait for it to be ready:

```bash
source continuous_integration/hdfs/startup_hdfs.sh
```

- Install dependencies and dask on the container

```bash
source continuous_integration/hdfs/install.sh
```

- Run the tests

```bash
source continuous_integration/hdfs/run_tests.sh
```

- Alternatively, you can start a terminal on the container and run the tests
  manually. This can be nicer for debugging:

```bash
# CONTAINER_ID should be defined from above, but if it isn't you can get it from
export CONTAINER_ID=$(docker ps -l -q)

# Start the bash session
docker exec -it $CONTAINER_ID bash

# Test just the hdfs tests
py.test dask/bytes/tests/test_hdfs.py -s -vv
```
There are a variety of other projects related to dask that are often
co-released.  We may want to check their status while releasing


Release per project:

*   Raise an issue in the https://github.com/dask/community issue tracker
    signaling your intent to release and the motivation.  Let that issue
    collect comments for a day to ensure that other maintainers are comfortable
    with releasing.

*   Update release notes in docs/source/changelog.rst

*   Commit

        git commit -a -m "bump version to x.x.x"

*   Tag commit

        git tag -a x.x.x -m 'Version x.x.x'

*   Push to GitHub

        git push dask main --tags

*   Upload to PyPI

        git clean -xfd
        python setup.py sdist bdist_wheel
        twine upload dist/*

*   Wait for [conda-forge](https://conda-forge.github.io) bots to track the
    change to PyPI

    This will typically happen in an hour or two.  There will be two PRs, one
    to `dask-core`, which you will likely be able to merge after tests pass,
    and another to `dask`, for which you might have to change version numbers
    if other packages (like `distributed`) are changing at the same time.

    In some cases you may also have to zero out build numbers.

    If for some reason you have to do this manually, then follow these steps:

    *  Update conda-smithy and run conda-smithy rerender

            git clone git@github.com:conda-forge/dask-core-feedstock
            cd dask-core-feedstock
            conda install conda-smithy
            conda-smithy rerender

    *  Get sha256 hash from pypi.org
    *  Update version number and hash in recipe
    *  Check dependencies
    *  Do the same for the dask-feedstock meta-package

*   Automated systems internal to Anaconda Inc then handle updating the
    Anaconda defaults channel
Dask
====

|Build Status| |Coverage| |Doc Status| |Discourse| |Version Status| |NumFOCUS|

Dask is a flexible parallel computing library for analytics.  See
documentation_ for more information.


LICENSE
-------

New BSD. See `License File <https://github.com/dask/dask/blob/main/LICENSE.txt>`__.

.. _documentation: https://dask.org
.. |Build Status| image:: https://github.com/dask/dask/workflows/CI/badge.svg?branch=main
   :target: https://github.com/dask/dask/actions?query=workflow%3A%22CI%22
.. |Coverage| image:: https://codecov.io/gh/dask/dask/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/dask/dask/branch/main
   :alt: Coverage status
.. |Doc Status| image:: https://readthedocs.org/projects/dask/badge/?version=latest
   :target: https://dask.org
   :alt: Documentation Status
.. |Discourse| image:: https://img.shields.io/discourse/users?logo=discourse&server=https%3A%2F%2Fdask.discourse.group
   :alt: Discuss Dask-related things and ask for help
   :target: https://dask.discourse.group
.. |Version Status| image:: https://img.shields.io/pypi/v/dask.svg
   :target: https://pypi.python.org/pypi/dask/
.. |NumFOCUS| image:: https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A
   :target: https://www.numfocus.org/
To build a local copy of the Dask documentation, install the packages in
``requirements-docs.txt`` and run ``make html``.

Optionally create and activate a ``conda`` environment first::

  conda create -n daskdocs -c conda-forge python=3.8
  conda activate daskdocs

Install the dependencies with ``pip``::

  python -m pip install -r requirements-docs.txt

After running ``make html`` the generated HTML documentation can be found in
the ``build/html`` directory. Open ``build/html/index.html`` to view the home
page for the documentation.
API
===

.. currentmodule:: dask.bag

Create Bags
-----------

.. autosummary::
   :toctree: generated/

   from_sequence
   from_delayed
   from_url
   range
   read_text
   read_avro

From dataframe
~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe

.. autosummary::
   :toctree: generated/

   DataFrame.to_bag
   Series.to_bag

Top-level functions
-------------------

.. currentmodule:: dask.bag

.. autosummary::
   :toctree: generated/

   concat
   map
   map_partitions
   to_textfiles
   zip

Random Sampling
---------------

.. autosummary::
   :toctree: generated/

   random.choices
   random.sample


Turn Bags into other things
---------------------------

.. autosummary::
   :toctree: generated/

   Bag.to_textfiles
   Bag.to_dataframe
   Bag.to_delayed
   Bag.to_avro


Bag Methods
-----------

.. autosummary::
   :toctree: generated/

   Bag
   Bag.accumulate
   Bag.all
   Bag.any
   Bag.compute
   Bag.count
   Bag.distinct
   Bag.filter
   Bag.flatten
   Bag.fold
   Bag.foldby
   Bag.frequencies
   Bag.groupby
   Bag.join
   Bag.map
   Bag.map_partitions
   Bag.max
   Bag.mean
   Bag.min
   Bag.persist
   Bag.pluck
   Bag.product
   Bag.reduction
   Bag.random_sample
   Bag.remove
   Bag.repartition
   Bag.starmap
   Bag.std
   Bag.sum
   Bag.take
   Bag.to_avro
   Bag.to_dataframe
   Bag.to_delayed
   Bag.to_textfiles
   Bag.topk
   Bag.var
   Bag.visualize


Item Methods
------------

.. autosummary::
   :toctree: generated/

   Item
   Item.apply
   Item.compute
   Item.from_delayed
   Item.persist
   Item.to_delayed
   Item.visualize
Bag
===

.. toctree::
   :maxdepth: 1
   :hidden:

   bag-creation.rst

Dask Bag implements operations like ``map``, ``filter``, ``fold``, and
``groupby`` on collections of generic Python objects.  It does this in parallel with a
small memory footprint using Python iterators.  It is similar to a parallel
version of PyToolz_ or a Pythonic version of the `PySpark RDD`_.

.. _PyToolz: https://toolz.readthedocs.io/en/latest/
.. _`PySpark RDD`: https://spark.apache.org/docs/latest/api/python

Examples
--------

Visit https://examples.dask.org/bag.html to see and run examples using Dask Bag.

Design
------

Dask bags coordinate many Python lists or Iterators, each of which forms a
partition of a larger collection.

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/-qIiJ1XtSv0"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

Common Uses
-----------

Dask bags are often used to parallelize simple computations on unstructured or
semi-structured data like text data, log files, JSON records, or user defined
Python objects.

Execution
---------

Execution on bags provide two benefits:

1.  Parallel: data is split up, allowing multiple cores or machines to execute
    in parallel
2.  Iterating: data processes lazily, allowing smooth execution of
    larger-than-memory data, even on a single machine within a single partition


Default scheduler
~~~~~~~~~~~~~~~~~

By default, ``dask.bag`` uses ``dask.multiprocessing`` for computation.  As a
benefit, Dask bypasses the GIL_ and uses multiple cores on pure Python objects.
As a drawback, Dask Bag doesn't perform well on computations that include a
great deal of inter-worker communication.  For common operations this is rarely
an issue as most Dask Bag workflows are embarrassingly parallel or result in
reductions with little data moving between workers.

.. _GIL: https://docs.python.org/3/glossary.html#term-gil


Shuffle
~~~~~~~

Some operations, like ``groupby``, require substantial inter-worker
communication. On a single machine, Dask uses partd_ to perform efficient,
parallel, spill-to-disk shuffles. When working in a cluster, Dask uses a task
based shuffle.

These shuffle operations are expensive and better handled by projects like
``dask.dataframe``. It is best to use ``dask.bag`` to clean and process data,
then transform it into an array or DataFrame before embarking on the more
complex operations that require shuffle steps.

.. _partd: https://github.com/mrocklin/partd


Known Limitations
-----------------

Bags provide very general computation (any Python function).  This generality
comes at cost.  Bags have the following known limitations:

1.  By default, they rely on the multiprocessing scheduler, which has its own
    set of known limitations (see :doc:`shared`)
2.  Bags are immutable and so you can not change individual elements
3.  Bag operations tend to be slower than array/DataFrame computations in the
    same way that standard Python containers tend to be slower than NumPy
    arrays and Pandas DataFrames
4.  Bag's ``groupby`` is slow.  You should try to use Bag's ``foldby`` if possible.
    Using ``foldby`` requires more thought though


Name
----

*Bag* is the mathematical name for an unordered collection allowing repeats. It
is a friendly synonym to multiset_. A bag, or a multiset, is a generalization of
the concept of a set that, unlike a set, allows multiple instances of the
multiset's elements:

* ``list``: *ordered* collection *with repeats*, ``[1, 2, 3, 2]``
* ``set``: *unordered* collection *without repeats*,  ``{1, 2, 3}``
* ``bag``: *unordered* collection *with repeats*, ``{1, 2, 2, 3}``

So, a bag is like a list, but it doesn't guarantee an ordering among elements.
There can be repeated elements but you can't ask for the ith element.

.. _multiset: https://en.wikipedia.org/wiki/Bag_(mathematics)
Custom Collections
==================

For many problems, the built-in Dask collections (``dask.array``,
``dask.dataframe``, ``dask.bag``, and ``dask.delayed``) are sufficient. For
cases where they aren't, it's possible to create your own Dask collections. Here
we describe the required methods to fulfill the Dask collection interface.

.. note:: This is considered an advanced feature. For most cases the built-in
          collections are probably sufficient.

Before reading this you should read and understand:

- :doc:`overview <graphs>`
- :doc:`graph specification <spec>`
- :doc:`custom graphs <custom-graphs>`

**Contents**

- :ref:`Description of the Dask collection interface <collection-interface>`
- :ref:`How this interface is used to implement the core Dask
  methods <core-method-internals>`
- :ref:`How to add the core methods to your class <adding-methods-to-class>`
- :ref:`example-dask-collection`
- :ref:`How to check if something is a Dask collection <is-dask-collection>`
- :ref:`How to make tokenize work with your collection <deterministic-hashing>`


.. _collection-interface:

The Dask Collection Interface
-----------------------------

To create your own Dask collection, you need to fulfill the following
interface. Note that there is no required base class.

It is recommended to also read :ref:`core-method-internals` to see how this
interface is used inside Dask.


.. method:: __dask_graph__(self)

    The Dask graph.

    **Returns**

    dsk : Mapping, None
        The Dask graph.  If ``None``, this instance will not be interpreted as a
        Dask collection, and none of the remaining interface methods will be
        called.

    If the collection also specifies :meth:`__dask_layers__`, then ``dsk`` must be a
    :class:`~dask.highlevelgraph.HighLevelGraph` or ``None``.


.. method:: __dask_keys__(self)

    The output keys for the Dask graph.

    Note that there are additional constraints on keys for a Dask collection
    than those described in the :doc:`task graph specification documentation <spec>`.
    These additional constraints are described below.

    **Returns**

    keys : list
        A possibly nested list of keys that represent the outputs of the graph.
        After computation, the results will be returned in the same layout,
        with the keys replaced with their corresponding outputs.

    All keys must either be non-empty strings or tuples where the first element is a
    non-empty string, followed by zero or more arbitrary hashables.
    The non-empty string is commonly known as the *collection name*. All collections
    embedded in the dask package have exactly one name, but this is not a requirement.

    These are all valid outputs:

    - ``[]``
    - ``["x", "y"]``
    - ``[[("y", "a", 0), ("y", "a", 1)], [("y", "b", 0), ("y", "b", 1)]``


.. method:: __dask_layers__(self)

    This method only needs to be implemented if the collection uses
    :class:`~dask.highlevelgraph.HighLevelGraph` to implement its dask graph.

    **Returns**

    names : tuple
        Tuple of names of the HighLevelGraph layers which contain all keys returned by
        :meth:`__dask_keys__`.


.. staticmethod:: __dask_optimize__(dsk, keys, \*\*kwargs)

    Given a graph and keys, return a new optimized graph.

    This method can be either a ``staticmethod`` or a ``classmethod``, but not
    an ``instancemethod``.

    Note that graphs and keys are merged before calling ``__dask_optimize__``;
    as such, the graph and keys passed to this method may represent more than
    one collection sharing the same optimize method.

    If not implemented, defaults to returning the graph unchanged.

    **Parameters**

    dsk : MutableMapping
        The merged graphs from all collections sharing the same
        ``__dask_optimize__`` method.
    keys : list
        A list of the outputs from ``__dask_keys__`` from all collections
        sharing the same ``__dask_optimize__`` method.
    \*\*kwargs
        Extra keyword arguments forwarded from the call to ``compute`` or
        ``persist``. Can be used or ignored as needed.

    **Returns**

    optimized_dsk : MutableMapping
        The optimized Dask graph.


.. staticmethod:: __dask_scheduler__(dsk, keys, \*\*kwargs)

    The default scheduler ``get`` to use for this object.

    Usually attached to the class as a staticmethod, e.g.:

    >>> import dask.threaded
    >>> class MyCollection:
    ...     # Use the threaded scheduler by default
    ...     __dask_scheduler__ = staticmethod(dask.threaded.get)


.. method:: __dask_postcompute__(self)

    Return the finalizer and (optional) extra arguments to convert the computed
    results into their in-memory representation.

    Used to implement ``dask.compute``.

    **Returns**

    finalize : callable
        A function with the signature ``finalize(results, *extra_args)``.
        Called with the computed results in the same structure as the
        corresponding keys from ``__dask_keys__``, as well as any extra
        arguments as specified in ``extra_args``. Should perform any necessary
        finalization before returning the corresponding in-memory collection
        from ``compute``. For example, the ``finalize`` function for
        ``dask.array.Array`` concatenates all the individual array chunks into
        one large numpy array, which is then the result of ``compute``.
    extra_args : tuple
        Any extra arguments to pass to ``finalize`` after ``results``. If no
        extra arguments should be an empty tuple.


.. method:: __dask_postpersist__(self)

    Return the rebuilder and (optional) extra arguments to rebuild an equivalent
    Dask collection from a persisted or rebuilt graph.

    Used to implement :func:`dask.persist`.

    **Returns**

    rebuild : callable
        A function with the signature
        ``rebuild(dsk, *extra_args, rename : Mapping[str, str] = None)``.
        ``dsk`` is a Mapping which contains at least the output keys returned by
        :meth:`__dask_keys__`. The callable should return an equivalent Dask collection
        with the same keys as ``self``, but with the results that are computed through a
        different graph. In the case of :func:`dask.persist`, the new graph will have
        just the output keys and the values already computed.

        If the optional parameter ``rename`` is specified, it indicates that output
        keys may be changing too; e.g. if the previous output of :meth:`__dask_keys__`
        was ``[("a", 0), ("a", 1)]``, after calling
        ``rebuild(dsk, *extra_args, rename={"a": "b"})`` it must become
        ``[("b", 0), ("b", 1)]``.
        The ``rename`` mapping may not contain the collection name(s); in such case the
        associated keys do not change. It may contain replacements for unexpected names,
        which must be ignored.

    extra_args : tuple
        Any extra arguments to pass to ``rebuild`` after ``dsk``. If no extra
        arguments are necessary, it must be an empty tuple.


.. note:: It's also recommended to define ``__dask_tokenize__``,
          see :ref:`deterministic-hashing`.


.. _core-method-internals:

Internals of the Core Dask Methods
----------------------------------

Dask has a few *core* functions (and corresponding methods) that implement
common operations:

- ``compute``: Convert one or more Dask collections into their in-memory
  counterparts
- ``persist``: Convert one or more Dask collections into equivalent Dask
  collections with their results already computed and cached in memory
- ``optimize``: Convert one or more Dask collections into equivalent Dask
  collections sharing one large optimized graph
- ``visualize``: Given one or more Dask collections, draw out the graph that
  would be passed to the scheduler during a call to ``compute`` or ``persist``

Here we briefly describe the internals of these functions to illustrate how
they relate to the above interface.

Compute
~~~~~~~

The operation of ``compute`` can be broken into three stages:

1. **Graph Merging & Optimization**

   First, the individual collections are converted to a single large graph and
   nested list of keys. How this happens depends on the value of the
   ``optimize_graph`` keyword, which each function takes:

   - If ``optimize_graph`` is ``True`` (default), then the collections are first
     grouped by their ``__dask_optimize__`` methods.  All collections with the
     same ``__dask_optimize__`` method have their graphs merged and keys
     concatenated, and then a single call to each respective
     ``__dask_optimize__`` is made with the merged graphs and keys.  The
     resulting graphs are then merged.

   - If ``optimize_graph`` is ``False``, then all the graphs are merged and all
     the keys concatenated.

   After this stage there is a single large graph and nested list of keys which
   represents all the collections.

2. **Computation**

   After the graphs are merged and any optimizations performed, the resulting
   large graph and nested list of keys are passed on to the scheduler.  The
   scheduler to use is chosen as follows:

   - If a ``get`` function is specified directly as a keyword, use that
   - Otherwise, if a global scheduler is set, use that
   - Otherwise fall back to the default scheduler for the given collections.
     Note that if all collections don't share the same ``__dask_scheduler__``
     then an error will be raised.

   Once the appropriate scheduler ``get`` function is determined, it is called
   with the merged graph, keys, and extra keyword arguments.  After this stage,
   ``results`` is a nested list of values. The structure of this list mirrors
   that of ``keys``, with each key substituted with its corresponding result.

3. **Postcompute**

   After the results are generated, the output values of ``compute`` need to be
   built. This is what the ``__dask_postcompute__`` method is for.
   ``__dask_postcompute__`` returns two things:

   - A ``finalize`` function, which takes in the results for the corresponding
     keys
   - A tuple of extra arguments to pass to ``finalize`` after the results

   To build the outputs, the list of collections and results is iterated over,
   and the finalizer for each collection is called on its respective results.

In pseudocode, this process looks like the following:

.. code:: python

    def compute(*collections, **kwargs):
        # 1. Graph Merging & Optimization
        # -------------------------------
        if kwargs.pop('optimize_graph', True):
            # If optimization is turned on, group the collections by
            # optimization method, and apply each method only once to the merged
            # sub-graphs.
            optimization_groups = groupby_optimization_methods(collections)
            graphs = []
            for optimize_method, cols in optimization_groups:
                # Merge the graphs and keys for the subset of collections that
                # share this optimization method
                sub_graph = merge_graphs([x.__dask_graph__() for x in cols])
                sub_keys = [x.__dask_keys__() for x in cols]
                # kwargs are forwarded to ``__dask_optimize__`` from compute
                optimized_graph = optimize_method(sub_graph, sub_keys, **kwargs)
                graphs.append(optimized_graph)
            graph = merge_graphs(graphs)
        else:
            graph = merge_graphs([x.__dask_graph__() for x in collections])
        # Keys are always the same
        keys = [x.__dask_keys__() for x in collections]

        # 2. Computation
        # --------------
        # Determine appropriate get function based on collections, global
        # settings, and keyword arguments
        get = determine_get_function(collections, **kwargs)
        # Pass the merged graph, keys, and kwargs to ``get``
        results = get(graph, keys, **kwargs)

        # 3. Postcompute
        # --------------
        output = []
        # Iterate over the results and collections
        for res, collection in zip(results, collections):
            finalize, extra_args = collection.__dask_postcompute__()
            out = finalize(res, **extra_args)
            output.append(out)

        # `dask.compute` always returns tuples
        return tuple(output)


Persist
~~~~~~~

Persist is very similar to ``compute``, except for how the return values are
created. It too has three stages:

1. **Graph Merging & Optimization**

   Same as in ``compute``.

2. **Computation**

   Same as in ``compute``, except in the case of the distributed scheduler,
   where the values in ``results`` are futures instead of values.

3. **Postpersist**

   Similar to ``__dask_postcompute__``, ``__dask_postpersist__`` is used to
   rebuild values in a call to ``persist``. ``__dask_postpersist__`` returns
   two things:

   - A ``rebuild`` function, which takes in a persisted graph.  The keys of
     this graph are the same as ``__dask_keys__`` for the corresponding
     collection, and the values are computed results (for the single-machine
     scheduler) or futures (for the distributed scheduler).
   - A tuple of extra arguments to pass to ``rebuild`` after the graph

   To build the outputs of ``persist``, the list of collections and results is
   iterated over, and the rebuilder for each collection is called on the graph
   for its respective results.

In pseudocode, this looks like the following:

.. code:: python

    def persist(*collections, **kwargs):
        # 1. Graph Merging & Optimization
        # -------------------------------
        # **Same as in compute**
        graph = ...
        keys = ...

        # 2. Computation
        # --------------
        # **Same as in compute**
        results = ...

        # 3. Postpersist
        # --------------
        output = []
        # Iterate over the results and collections
        for res, collection in zip(results, collections):
            # res has the same structure as keys
            keys = collection.__dask_keys__()
            # Get the computed graph for this collection.
            # Here flatten converts a nested list into a single list
            subgraph = {k: r for (k, r) in zip(flatten(keys), flatten(res))}

            # Rebuild the output dask collection with the computed graph
            rebuild, extra_args = collection.__dask_postpersist__()
            out = rebuild(subgraph, *extra_args)

            output.append(out)

        # dask.persist always returns tuples
        return tuple(output)


Optimize
~~~~~~~~

The operation of ``optimize`` can be broken into two stages:

1. **Graph Merging & Optimization**

   Same as in ``compute``.

2. **Rebuilding**

   Similar to ``persist``, the ``rebuild`` function and arguments from
   ``__dask_postpersist__`` are used to reconstruct equivalent collections from
   the optimized graph.

In pseudocode, this looks like the following:

.. code:: python

    def optimize(*collections, **kwargs):
        # 1. Graph Merging & Optimization
        # -------------------------------
        # **Same as in compute**
        graph = ...

        # 2. Rebuilding
        # -------------
        # Rebuild each dask collection using the same large optimized graph
        output = []
        for collection in collections:
            rebuild, extra_args = collection.__dask_postpersist__()
            out = rebuild(graph, *extra_args)
            output.append(out)

        # dask.optimize always returns tuples
        return tuple(output)


Visualize
~~~~~~~~~

Visualize is the simplest of the 4 core functions. It only has two stages:

1. **Graph Merging & Optimization**

   Same as in ``compute``.

2. **Graph Drawing**

   The resulting merged graph is drawn using ``graphviz`` and outputs to the
   specified file.

In pseudocode, this looks like the following:

.. code:: python

    def visualize(*collections, **kwargs):
        # 1. Graph Merging & Optimization
        # -------------------------------
        # **Same as in compute**
        graph = ...

        # 2. Graph Drawing
        # ----------------
        # Draw the graph with graphviz's `dot` tool and return the result.
        return dot_graph(graph, **kwargs)


.. _adding-methods-to-class:

Adding the Core Dask Methods to Your Class
------------------------------------------

Defining the above interface will allow your object to used by the core Dask
functions (``dask.compute``, ``dask.persist``, ``dask.visualize``, etc.). To
add corresponding method versions of these, you can subclass from
``dask.base.DaskMethodsMixin`` which adds implementations of ``compute``,
``persist``, and ``visualize`` based on the interface above.

.. _example-dask-collection:

Example Dask Collection
-----------------------

Here we create a Dask collection representing a tuple.  Every element in the
tuple is represented as a task in the graph.  Note that this is for illustration
purposes only - the same user experience could be done using normal tuples with
elements of ``dask.delayed``:

.. code:: python

    # Saved as dask_tuple.py
    import dask
    from dask.base import DaskMethodsMixin, replace_name_in_key
    from dask.optimization import cull

    # We subclass from DaskMethodsMixin to add common dask methods to our
    # class. This is nice but not necessary for creating a dask collection.
    class Tuple(DaskMethodsMixin):
        def __init__(self, dsk, keys):
            # The init method takes in a dask graph and a set of keys to use
            # as outputs.
            self._dsk = dsk
            self._keys = keys

        def __dask_graph__(self):
            return self._dsk

        def __dask_keys__(self):
            return self._keys

        @staticmethod
        def __dask_optimize__(dsk, keys, **kwargs):
            # We cull unnecessary tasks here. Note that this isn't necessary,
            # dask will do this automatically, this just shows one optimization
            # you could do.
            dsk2, _ = cull(dsk, keys)
            return dsk2

        # Use the threaded scheduler by default.
        __dask_scheduler__ = staticmethod(dask.threaded.get)

        def __dask_postcompute__(self):
            # We want to return the results as a tuple, so our finalize
            # function is `tuple`. There are no extra arguments, so we also
            # return an empty tuple.
            return tuple, ()

        def __dask_postpersist__(self):
            # We need to return a callable with the signature
            # rebuild(dsk, *extra_args, rename: Mapping[str, str] = None)
            return Tuple._rebuild, (self._keys,)

        @staticmethod
        def _rebuild(dsk, keys, *, rename=None):
            if rename is not None:
                keys = [replace_name_in_key(key, rename) for key in keys]
            return Tuple(dsk, keys)

        def __dask_tokenize__(self):
            # For tokenize to work we want to return a value that fully
            # represents this object. In this case it's the list of keys
            # to be computed.
            return self._keys

Demonstrating this class:

.. code:: python

    >>> from dask_tuple import Tuple
    >>> from operator import add, mul

    # Define a dask graph
    >>> dsk = {"k0": 1,
    ...        ("x", "k1"): 2,
    ...        ("x", 1): (add, "k0", ("x", "k1")),
    ...        ("x", 2): (mul, ("x", "k1"), 2),
    ...        ("x", 3): (add, ("x", "k1"), ("x", 1))}

    # The output keys for this graph.
    # The first element of each tuple must be the same across the whole collection;
    # the remainder are arbitrary, unique hashables
    >>> keys = [("x", "k1"), ("x", 1), ("x", 2), ("x", 3)]

    >>> x = Tuple(dsk, keys)

    # Compute turns Tuple into a tuple
    >>> x.compute()
    (2, 3, 4, 5)

    # Persist turns Tuple into a Tuple, with each task already computed
    >>> x2 = x.persist()
    >>> isinstance(x2, Tuple)
    True
    >>> x2.__dask_graph__()
    {('x', 'k1'): 2, ('x', 1): 3, ('x', 2): 4, ('x', 3): 5}
    >>> x2.compute()
    (2, 3, 4, 5)


.. _is-dask-collection:

Checking if an object is a Dask collection
------------------------------------------

To check if an object is a Dask collection, use
``dask.base.is_dask_collection``:

.. code:: python

    >>> from dask.base import is_dask_collection
    >>> from dask import delayed

    >>> x = delayed(sum)([1, 2, 3])
    >>> is_dask_collection(x)
    True
    >>> is_dask_collection(1)
    False


.. _deterministic-hashing:

Implementing Deterministic Hashing
----------------------------------

Dask implements its own deterministic hash function to generate keys based on
the value of arguments.  This function is available as ``dask.base.tokenize``.
Many common types already have implementations of ``tokenize``, which can be
found in ``dask/base.py``.

When creating your own custom classes, you may need to register a ``tokenize``
implementation. There are two ways to do this:

1. The ``__dask_tokenize__`` method

   Where possible, it is recommended to define the ``__dask_tokenize__`` method.
   This method takes no arguments and should return a value fully
   representative of the object.

2. Register a function with ``dask.base.normalize_token``

   If defining a method on the class isn't possible or you need to
   customize the tokenize function for a class whose super-class is
   already registered (for example if you need to sub-class built-ins),
   you can register a tokenize function with the ``normalize_token``
   dispatch.  The function should have the same signature as described
   above.

In both cases the implementation should be the same, where only the location of the
definition is different.

.. note:: Both Dask collections and normal Python objects can have
          implementations of ``tokenize`` using either of the methods
          described above.

Example
~~~~~~~

.. code:: python

    >>> from dask.base import tokenize, normalize_token

    # Define a tokenize implementation using a method.
    >>> class Foo:
    ...     def __init__(self, a, b):
    ...         self.a = a
    ...         self.b = b
    ...
    ...     def __dask_tokenize__(self):
    ...         # This tuple fully represents self
    ...         return (Foo, self.a, self.b)

    >>> x = Foo(1, 2)
    >>> tokenize(x)
    '5988362b6e07087db2bc8e7c1c8cc560'
    >>> tokenize(x) == tokenize(x)  # token is deterministic
    True

    # Register an implementation with normalize_token
    >>> class Bar:
    ...     def __init__(self, x, y):
    ...         self.x = x
    ...         self.y = y

    >>> @normalize_token.register(Bar)
    ... def tokenize_bar(x):
    ...     return (Bar, x.x, x.x)

    >>> y = Bar(1, 2)
    >>> tokenize(y)
    '5a7e9c3645aa44cf13d021c14452152e'
    >>> tokenize(y) == tokenize(y)
    True
    >>> tokenize(y) == tokenize(x)  # tokens for different objects aren't equal
    False


For more examples, see ``dask/base.py`` or any of the built-in Dask collections.
.. _dataframe.performance:

Best Practices
==============

It is easy to get started with Dask DataFrame, but using it *well* does require
some experience.  This page contains suggestions for best practices, and
includes solutions to common problems.


Use Pandas
----------

For data that fits into RAM, Pandas can often be faster and easier to use than
Dask DataFrame.  While "Big Data" tools can be exciting, they are almost always
worse than normal data tools while those remain appropriate.


Reduce, and then use Pandas
---------------------------

Similar to above, even if you have a large dataset there may be a point in your
computation where you've reduced things to a more manageable level.  You may
want to switch to Pandas at this point.

.. code-block:: python

   df = dd.read_parquet('my-giant-file.parquet')
   df = df[df.name == 'Alice']              # Select a subsection
   result = df.groupby('id').value.mean()   # Reduce to a smaller size
   result = result.compute()                # Convert to Pandas dataframe
   result...                                # Continue working with Pandas

Pandas Performance Tips Apply to Dask DataFrame
-----------------------------------------------

Usual Pandas performance tips like avoiding apply, using vectorized
operations, using categoricals, etc., all apply equally to Dask DataFrame.  See
`Modern Pandas <https://tomaugspurger.github.io/modern-1-intro>`_ by `Tom
Augspurger <https://github.com/TomAugspurger>`_ for a good read on this topic.

Use the Index
-------------

Dask DataFrame can be optionally sorted along a single index column.  Some
operations against this column can be very fast.  For example, if your dataset
is sorted by time, you can quickly select data for a particular day, perform
time series joins, etc.  You can check if your data is sorted by looking at the
``df.known_divisions`` attribute.  You can set an index column using the
``.set_index(column_name)`` method.  This operation is expensive though, so use
it sparingly (see below):

.. code-block:: python

   df = df.set_index('timestamp')  # set the index to make some operations fast

   df.loc['2001-01-05':'2001-01-12']  # this is very fast if you have an index
   df.merge(df2, left_index=True, right_index=True)  # this is also very fast

For more information, see documentation on :ref:`dataframe partitions <dataframe-design-partitions>`.


Avoid Full-Data Shuffling
-------------------------

Setting an index is an important but expensive operation (see above).  You
should do it infrequently and you should persist afterwards (see below).

Some operations like ``set_index`` and ``merge/join`` are harder to do in a
parallel or distributed setting than if they are in-memory on a single machine.
In particular, *shuffling operations* that rearrange data become much more
communication intensive.  For example, if your data is arranged by customer ID
but now you want to arrange it by time, all of your partitions will have to talk
to each other to exchange shards of data.  This can be an intensive process,
particularly on a cluster.

So, definitely set the index but try do so infrequently.  After you set the
index, you may want to ``persist`` your data if you are on a cluster:

.. code-block:: python

   df = df.set_index('column_name')  # do this infrequently

Additionally, ``set_index`` has a few options that can accelerate it in some
situations.  For example, if you know that your dataset is sorted or you already
know the values by which it is divided, you can provide these to accelerate the
``set_index`` operation.  For more information, see the `set_index docstring
<https://docs.dask.org/en/latest/dataframe-api.html#dask.dataframe.DataFrame.set_index>`_.

.. code-block:: python

   df2 = df.set_index(d.timestamp, sorted=True)


Persist Intelligently
---------------------

.. note:: This section is only relevant to users on distributed systems.

Often DataFrame workloads look like the following:

1.  Load data from files
2.  Filter data to a particular subset
3.  Shuffle data to set an intelligent index
4.  Several complex queries on top of this indexed data

It is often ideal to load, filter, and shuffle data once and keep this result in
memory.  Afterwards, each of the several complex queries can be based off of
this in-memory data rather than have to repeat the full load-filter-shuffle
process each time.  To do this, use the `client.persist
<https://distributed.dask.org/en/latest/api.html#distributed.Client.persist>`_
method:

.. code-block:: python

   df = dd.read_csv('s3://bucket/path/to/*.csv')
   df = df[df.balance < 0]
   df = client.persist(df)

   df = df.set_index('timestamp')
   df = client.persist(df)

   >>> df.customer_id.nunique().compute()
   18452844

   >>> df.groupby(df.city).size().compute()
   ...

Persist is important because Dask DataFrame is *lazy by default*.  It is a
way of telling the cluster that it should start executing the computations
that you have defined so far, and that it should try to keep those results in
memory.  You will get back a new DataFrame that is semantically equivalent to
your old DataFrame, but now points to running data.  Your old DataFrame still
points to lazy computations:

.. code-block:: python

   # Don't do this
   client.persist(df)  # persist doesn't change the input in-place

   # Do this instead
   df = client.persist(df)  # replace your old lazy DataFrame


Repartition to Reduce Overhead
------------------------------

Your Dask DataFrame is split up into many Pandas DataFrames.  We sometimes call
these "partitions", and often the number of partitions is decided for you. For
example, it might be the number of CSV files from which you are reading. However,
over time, as you reduce or increase the size of your pandas DataFrames by
filtering or joining, it may be wise to reconsider how many partitions you need.
There is a cost to having too many or having too few.

Partitions should fit comfortably in memory (smaller than a gigabyte) but also
not be too many.  Every operation on every partition takes the central
scheduler a few hundred microseconds to process.  If you have a few thousand
tasks this is barely noticeable, but it is nice to reduce the number if
possible.

A common situation is that you load lots of data into reasonably sized
partitions (Dask's defaults make decent choices), but then you filter down your
dataset to only a small fraction of the original.  At this point, it is wise to
regroup your many small partitions into a few larger ones.  You can do this by
using the :py:class:`dask.dataframe.DataFrame.repartition` method:

.. code-block:: python

   df = dd.read_csv('s3://bucket/path/to/*.csv')
   df = df[df.name == 'Alice']  # only 1/100th of the data
   df = df.repartition(npartitions=df.npartitions // 100)

   df = df.persist()  # if on a distributed system

This helps to reduce overhead and increase the effectiveness of vectorized
Pandas operations.  You should aim for partitions that have around 100MB of
data each.

Additionally, reducing partitions is very helpful just before shuffling, which
creates ``n log(n)`` tasks relative to the number of partitions.  DataFrames
with less than 100 partitions are much easier to shuffle than DataFrames with
tens of thousands.


Joins
-----

Joining two DataFrames can be either very expensive or very cheap depending on
the situation.  It is cheap in the following cases:

1.  Joining a Dask DataFrame with a Pandas DataFrame
2.  Joining a Dask DataFrame with another Dask DataFrame of a single partition
3.  Joining Dask DataFrames along their indexes

Also, it is expensive in the following case:

1.  Joining Dask DataFrames along columns that are not their index

The expensive case requires a shuffle.  This is fine, and Dask DataFrame will
complete the job well, but it will be more expensive than a typical linear-time
operation:

.. code-block:: python

   dd.merge(a, pandas_df)  # fast
   dd.merge(a, b, left_index=True, right_index=True)  # fast
   dd.merge(a, b, left_index=True, right_on='id')  # half-fast, half-slow
   dd.merge(a, b, left_on='id', right_on='id')  # slow

For more information see :doc:`Joins <dataframe-joins>`.


Store Data in Apache Parquet Format
-----------------------------------

HDF5 is a popular choice for Pandas users with high performance needs.  We
encourage Dask DataFrame users to :doc:`store and load data <dataframe-create>`
using Parquet instead.  `Apache Parquet <https://parquet.apache.org/>`_ is a
columnar binary format that is easy to split into multiple files (easier for
parallel loading) and is generally much simpler to deal with than HDF5 (from
the library's perspective).  It is also a common format used by other big data
systems like `Apache Spark <https://spark.apache.org/>`_ and `Apache Impala
<https://impala.apache.org/>`_, and so it is useful to interchange with other
systems:

.. code-block:: python

   df.to_parquet('path/to/my-results/')
   df = dd.read_parquet('path/to/my-results/')

Dask supports reading parquet files with different engine implementations of
the Apache Parquet format for Python:

.. code-block:: python

   df1 = dd.read_parquet('path/to/my-results/', engine='fastparquet')
   df2 = dd.read_parquet('path/to/my-results/', engine='pyarrow')

These libraries can be installed using:

.. code-block:: shell

   conda install fastparquet pyarrow -c conda-forge

`fastparquet <https://github.com/dask/fastparquet/>`_ is a Python-based
implementation that uses the `Numba <https://numba.pydata.org/>`_
Python-to-LLVM compiler. PyArrow is part of the
`Apache Arrow <https://arrow.apache.org/>`_ project and uses the `C++
implementation of Apache Parquet <https://github.com/apache/parquet-cpp>`_.
:orphan:

Images and Logos
================

.. image:: images/dask_icon.svg
   :alt: Dask logo.

.. image:: images/dask_icon_no_pad.svg
   :alt: Dask logo without padding.

.. image:: images/dask_horizontal.svg
   :alt: Dask logo.

.. image:: images/dask_horizontal_no_pad.svg
   :alt: Dask logo without padding.

.. image:: images/dask_horizontal_white.svg
   :alt: Dask logo.

.. image:: images/dask_horizontal_white_no_pad.svg
   :alt: Dask logo.

.. image:: images/dask_stacked.svg
   :alt: Dask logo.

.. image:: images/dask_stacked_white.svg
   :alt: Dask logo.

.. image:: images/HHMI_Janelia_Color.png
   :alt: HHMI Janelia logo.
Shared Memory
=============

The asynchronous scheduler accepts any ``concurrent.futures.Executor``
instance. This includes instances of the ``ThreadPoolExecutor`` and
``ProcessPoolExecutor`` defined in the Python standard library as well as any
other subclass from a 3rd party library. Dask also defines its own
``SynchronousExecutor`` for that simply runs functions on the main thread
(useful for debugging).

Full dask ``get`` functions exist in each of ``dask.threaded.get``,
``dask.multiprocessing.get`` and ``dask.get`` respectively.


Policy
------

The asynchronous scheduler maintains indexed data structures that show which
tasks depend on which data, what data is available, and what data is waiting on
what tasks to complete before it can be released, and what tasks are currently
running.  It can update these in constant time relative to the number of total
and available tasks.  These indexed structures make the dask async scheduler
scalable to very many tasks on a single machine.

.. image:: images/async-embarrassing.gif
   :width: 50 %
   :align: right
   :alt: Embarrassingly parallel dask flow

To keep the memory footprint small, we choose to keep ready-to-run tasks in a
last-in-first-out stack such that the most recently made available tasks get
priority.  This encourages the completion of chains of related tasks before new
chains are started.  This can also be queried in constant time.


Performance
-----------

**tl;dr** The threaded scheduler overhead behaves roughly as follows:

*  200us overhead per task
*  10us startup time (if you wish to make a new ThreadPoolExecutor each time)
*  Constant scaling with number of tasks
*  Linear scaling with number of dependencies per task

Schedulers introduce overhead.  This overhead effectively limits the
granularity of our parallelism.  Below we measure overhead of the async
scheduler with different apply functions (threaded, sync, multiprocessing), and
under different kinds of load (embarrassingly parallel, dense communication).

The quickest/simplest test we can do it to use IPython's ``timeit`` magic:

.. ipython::

   In [1]: import dask.array as da

   In [2]: x = da.ones(1000, chunks=(2,)).sum()

   In [3]: len(x.dask)
   Out[3]: 1168

   In [4]: %timeit x.compute()
   92.1 ms ± 2.61 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)

So this takes ~90 microseconds per task.  About 100ms of this is from overhead:

.. ipython::

   In [5]: x = da.ones(1000, chunks=(1000,)).sum()
   
   In [6]: %timeit x.compute()
   1.18 ms ± 8.64 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)

There is some overhead from spinning up a ThreadPoolExecutor each time.
This may be mediated by using a global or contextual pool:

.. code-block:: python

   >>> from concurrent.futures import ThreadPoolExecutor
   >>> pool = ThreadPoolExecutor()
   >>> dask.config.set(pool=pool)  # set global ThreadPoolExecutor

   or

   >>> with dask.config.set(pool=pool)  # use ThreadPoolExecutor throughout with block
   ...     ...

We now measure scaling the number of tasks and scaling the density of the
graph:

.. image:: images/trivial.svg
   :width: 30 %
   :align: right
   :alt: Adding nodes

Linear scaling with number of tasks
```````````````````````````````````

As we increase the number of tasks in a graph, we see that the scheduling
overhead grows linearly.  The asymptotic cost per task depends on the scheduler.
The schedulers that depend on some sort of asynchronous pool have costs of a few
milliseconds and the single threaded schedulers have costs of a few microseconds.

.. figure:: images/scaling-nodes.png
   :alt: Graph depicting how well Dask scales with the number of nodes in the task graph. Graph shows the duration in seconds on the y-axis versus number of edges per task on the x-axis. The time to schedule the entire graph is constant initially, followed by a linear increase after roughly 500 tasks for multiprocessing and threaded schedulers and 10 tasks for async and core schedulers. The inverse is true for the cost per task, with a linear cost decrease, followed by more or less constant cost.
   
   Scheduling overhead for the entire graph (left) vs. per task (right)

.. image:: images/crosstalk.svg
   :width: 40 %
   :align: right
   :alt: Adding edges

Linear scaling with number of edges
```````````````````````````````````

As we increase the number of edges per task, the scheduling overhead
again increases linearly.

Note: Neither the naive core scheduler nor the multiprocessing scheduler
are good at workflows with non-trivial cross-task
communication; they have been removed from the plot.

.. figure:: images/scaling-edges.png
   :alt: Graph depicting how well Dask scales with the number of edges in the task graph. Graph shows the duration in seconds on the y-axis versus number of edges per task on the x-axis. As the number of edges increases from 0 to 100, the time to schedule the entire graph using the threaded scheduler goes from 2 to 8 seconds whereas using the async scheduler goes from 0 to 3 seconds. The cost per edge decreases up until about 10 edges, after which the cost plateaus for both the threaded and async schedulers, with the async scheduler being consistenly faster.
   
   Scheduling overhead of the entire graph (left) vs. per edge (right)

`Download scheduling script`_


Known Limitations
-----------------

The shared memory scheduler has some notable limitations:

1.  It works on a single machine
2.  The threaded scheduler is limited by the GIL on Python code, so if your
    operations are pure python functions, you should not expect a multi-core
    speedup
3.  The multiprocessing scheduler must serialize functions between workers,
    which can fail
4.  The multiprocessing scheduler must serialize data between workers and the
    central process, which can be expensive
5.  The multiprocessing scheduler cannot transfer data directly between worker
    processes; all data routes through the main process.



.. _`Download scheduling script`: https://github.com/dask/dask/tree/main/docs/source/scripts/scheduling.py
:orphan:

.. this page is refernenced from dask.org. It might move to there some day.

Why Dask?
=========

This document gives high-level motivation on why people choose to adopt Dask.

.. contents:: :local:

Python's role in Data Science
-----------------------------

Python has grown to become the dominant language both in data analytics and
general programming:

.. image:: images/growth_of_languages.png
   :alt: Graph showing the growth of major programming languages based on Stack Overflow’s question views in World Bank high-income countries. A line graph with time from 2012 to 2018 on the x-axis and percent of overall question views each month on the y-axis. Python’s question views increase from about 4% to about 11% from 2012 to 2018, reaching the popularity of Java and JavaScript.
   :width: 75%

This is fueled both by computational libraries like Numpy, Pandas, and
Scikit-Learn and by a wealth of libraries for visualization, interactive
notebooks, collaboration, and so forth.

.. image:: images/growth_of_libraries.png
   :alt: Graph showing the growth of major python packages based on Stack Overflow's question views in World Bank high-income countries. A line graph with time on the x-axis from 2012 to 2018 and percent of overall question views each month on the y-axis. Pandas question views increased to about 0.9% in 2018, exceeding Django and NumPy.
   :width: 75%

However, these packages were not designed to scale beyond a single machine.
Dask was developed to scale these packages and the surrounding ecosystem.
It works with the existing Python ecosystem to scale it to multi-core
machines and distributed clusters.

*Image credit to Stack Overflow blogposts*
`#1 <https://stackoverflow.blog/2017/09/06/incredible-growth-python>`_
*and*
`#2 <https://stackoverflow.blog/2017/09/14/python-growing-quickly/>`_.


Dask has a Familiar API
-----------------------

Analysts often use tools like Pandas, Scikit-Learn, Numpy, and the rest of the
Python ecosystem to analyze data on their personal computer.  They like these
tools because they are efficient, intuitive, and widely trusted.  However, when
they choose to apply their analyses to larger datasets, they find that these
tools were not designed to scale beyond a single machine. And so, the analyst
rewrites their computation using a more scalable tool, often in another
language altogether.  This rewrite process slows down discovery and causes
frustration.

Dask provides ways to scale Pandas, Scikit-Learn, and Numpy workflows more
natively, with minimal rewriting.  It integrates well with these tools so that it copies
most of their API and uses their data structures internally.  Moreover, Dask is
co-developed with these libraries to ensure that they evolve consistently,
minimizing friction when transitioning from a local laptop,
to a multi-core workstation, and then to a distributed cluster.  Analysts
familiar with Pandas/Scikit-Learn/Numpy will be immediately familiar with their
Dask equivalents, and have much of their intuition carry over to a scalable
context.


Dask Scales out to Clusters
---------------------------

As datasets and computations scale faster than CPUs and RAM, we need to find
ways to scale our computations across multiple machines.  This introduces many
new concerns:

-  How to have computers talk to each other over the network?
-  How and when to move data between machines?
-  How to recover from machine failures?
-  How to deploy on an in-house cluster?
-  How to deploy on the cloud?
-  How to deploy on an HPC super-computer?
-  How to provide an API to this system that users find intuitive?
-  ...

While it is possible to build these systems in-house (and indeed, many exist),
many organizations increasingly depend on solutions developed within the
open source community.  These tend to be more robust, secure, and fully
featured without being tended by in-house staff.

Dask solves the problems above.  It figures out how to break up large
computations and route parts of them efficiently onto distributed hardware.
Dask is routinely run on thousand-machine clusters to process hundreds of
terabytes of data efficiently within secure environments.

Dask has utilities and documentation on how to deploy in-house, on
the cloud, or on HPC super-computers.  It supports encryption and
authentication using TLS/SSL certificates.  It is resilient and can handle the
failure of worker nodes gracefully and is elastic, and so can take advantage of
new nodes added on-the-fly.  Dask includes several user APIs that are used and
smoothed over by thousands of researchers across the globe working in different
domains.


Dask Scales Down to Single Computers
------------------------------------

*But a massive cluster is not always the right choice*

Today's laptops and workstations are surprisingly powerful and, if used
correctly, can handle datasets and computations for which we previously
depended on clusters.  A modern laptop has a multi-core CPU, 32GB of RAM, and
flash-based hard drives that can stream through data several times faster than
HDDs or SSDs of even a year or two ago.

As a result, Dask can empower analysts to manipulate 100GB+ datasets on their
laptop or 1TB+ datasets on a workstation without bothering with the cluster at
all.  This can be preferable for the following reasons:

1.  They can use their local software environment, rather than being
    constrained by what is available on the cluster or having to manage
    Docker images.
2.  They can more easily work while in transit, at a coffee shop, or at home
    away from the corporate network
3.  Debugging errors and analyzing performance is simpler and more pleasant on
    a single machine
4.  Their iteration cycles can be faster
5.  Their computations may be more efficient because all of the data is local
    and doesn't need to flow through the network or between separate processes

Dask can enable efficient parallel computations on single machines by
leveraging their multi-core CPUs and streaming data efficiently from disk.
It *can* run on a distributed cluster, but it doesn't *have* to.  Dask allows
you to swap out the cluster for single-machine schedulers which are surprisingly
lightweight, require no setup, and can run entirely within the same process as
the user's session.

To avoid excess memory use, Dask is good at finding ways to evaluate
computations in a low-memory footprint when possible by pulling in chunks of
data from disk, doing the necessary processing, and throwing away intermediate
values as quickly as possible.  This lets analysts perform computations on
moderately large datasets (100GB+) even on relatively low-power laptops.
This requires no configuration and no setup, meaning that adding Dask to a
single-machine computation adds very little cognitive overhead.

Dask is installed by default with `Anaconda <https://anaconda.com>`_
and so is already deployed on most data science machines.


Dask Integrates Natively with Python Code
-----------------------------------------

Python includes computational libraries like Numpy, Pandas, and Scikit-Learn,
and many others for data access, plotting, statistics, image and
signal processing, and more.  These libraries work together seamlessly to
produce a cohesive *ecosystem* of packages that co-evolve to meet the needs of
analysts in most domains today.

This ecosystem is tied together by common standards and protocols to which
everyone adheres, which allows these packages to benefit each other in
surprising and delightful ways.

Dask evolved from within this ecosystem.  It abides by these standards and
protocols and actively engages in community efforts to push forward new ones.
This enables the rest of the ecosystem to benefit from parallel and distributed
computing with minimal coordination.  Dask does not seek to disrupt or displace
the existing ecosystem, but rather to complement and benefit it from within.

As a result, Dask development is pushed forward by developer communities
from Pandas, Numpy, Scikit-Learn, Scikit-Image, Jupyter, and others.  This
engagement from the broader community growth helps users to trust the project
and helps to ensure that the Python ecosystem will continue to evolve in a
smooth and sustainable manner.


Dask Supports Complex Applications
----------------------------------

Some parallel computations are simple and just apply the same routine onto many
inputs without any kind of coordination.  These are simple to parallelize with
any system.

Somewhat more complex computations can be expressed with the
map-shuffle-reduce pattern popularized by Hadoop and Spark.
This is often sufficient to do most data cleaning tasks,
database-style queries, and some lightweight machine learning algorithms.

However, more complex parallel computations exist which do not fit into these
paradigms, and so are difficult to perform with traditional big-data
technologies.  These include more advanced algorithms for statistics or machine
learning, time series or local operations, or bespoke parallelism often found
within the systems of large enterprises.

Many companies and institutions today have problems which are
clearly parallelizable, but not clearly transformable into a big DataFrame
computation.  Today these companies tend to solve their problems either by
writing custom code with low-level systems like MPI, ZeroMQ, or sockets and
complex queuing systems, or by shoving their problem into a standard big-data
technology like MapReduce or Spark, and hoping for the best.

Dask helps to resolve these situations by exposing low-level APIs to its
internal task scheduler which is capable of executing very advanced
computations.  This gives engineers within the institution the ability to build
their own parallel computing system using the same engine that powers Dask's
arrays, DataFrames, and machine learning algorithms, but now with the
institution's own custom logic.  This allows engineers to keep complex
business logic in-house while still relying on Dask to handle network
communication, load balancing, resilience, diagnostics, etc..


Dask Delivers Responsive Feedback
---------------------------------

Because everything happens remotely, interactive parallel computing can be
frustrating for users.  They don't have a good sense of how computations are
progressing, what might be going wrong, or what parts of their code should they
focus on for performance.  The added distance between a user and their
computation can drastically affect how quickly they are able to identify and
resolve bugs and performance problems, which can drastically increase their
time to solution.

Dask keeps users informed and content with a suite of helpful diagnostic and
investigative tools including the following:

1.  A :doc:`real-time and responsive dashboard <understanding-performance>`
    that shows current progress, communication costs, memory use, and more,
    updated every 100ms
2.  A statistical profiler installed on every worker that polls each thread
    every 10ms to determine which lines in your code are taking up the most
    time across your entire computation
3.  An embedded IPython kernel in every worker and the scheduler, allowing
    users to directly investigate the state of their computation with a pop-up
    terminal
4.  The ability to reraise errors locally, so that they can use the traditional
    debugging tools to which they are accustomed, even when the error happens
    remotely


Links and More Information
--------------------------

From here you may want to read about some of our more common introductory
content:

-  :doc:`user-interfaces`
-  :doc:`scheduling`
-  :doc:`spark`
-  `Slides <https://dask.org/slides.html>`_
:orphan:

.. this page is referenced from the topbar which comes from the theme

Ecosystem
=========

There are a number of open source projects that extend the Dask interface and provide different
mechanisms for deploying Dask clusters. This is likely an incomplete list so if you spot something
missing - please `suggest a fix <https://github.com/dask/dask/edit/main/docs/source/ecosystem.rst>`_!

Building on Dask
----------------
Many packages include built-in support for Dask collections or wrap Dask collections internally
to enable parallelization.

Array
~~~~~
- `xarray <https://xarray.pydata.org>`_:  Wraps Dask
  Array, offering the same scalability, but with axis labels which add convenience when
  dealing with complex datasets.
- `cupy <https://docs.cupy.dev/en/stable>`_: Part of the Rapids project, GPU-enabled arrays
  can be used as the blocks of Dask Arrays. See the section :doc:`gpu` for more information.
- `sparse <https://github.com/pydata/sparse>`_: Implements sparse arrays of arbitrary dimension
  on top of ``numpy`` and ``scipy.sparse``.
- `pint <https://pint.readthedocs.io>`_: Allows arithmetic operations between them and conversions
  from and to different units.

DataFrame
~~~~~~~~~
- `cudf <https://docs.rapids.ai/api/cudf/stable/>`_: Part of the Rapids project, implements
  GPU-enabled dataframes which can be used as partitions in Dask Dataframes.
- `dask-geopandas <https://github.com/geopandas/dask-geopandas>`_: Early-stage subproject of
  geopandas, enabling parallelization of geopandas dataframes.

SQL
~~~
- `blazingSQL`_: Part of the Rapids project, implements SQL queries using ``cuDF``
  and Dask, for execution on CUDA/GPU-enabled hardware, including referencing
  externally-stored data.
- `dask-sql`_: Adds a SQL query layer on top of Dask.
  The API matches blazingSQL but it uses CPU instead of GPU. It still under development
  and not ready for a production use-case.
- `fugue-sql`_: Adds an abstract layer that makes code portable between across differing
  computing frameworks such as Pandas, Spark and Dask.

.. _blazingSQL: https://docs.blazingsql.com/
.. _dask-sql: https://dask-sql.readthedocs.io/en/latest/
.. _fugue-sql: https://fugue-tutorials.readthedocs.io/en/latest/tutorials/fugue_sql/index.html

Machine Learning
~~~~~~~~~~~~~~~~
- `dask-ml <https://ml.dask.org>`_: Implements distributed versions of common machine learning algorithms.
- `scikit-learn <https://scikit-learn.org/stable/>`_: Provide 'dask' to the joblib backend to parallelize
  scikit-learn algorithms with dask as the processor.
- `xgboost <https://xgboost.readthedocs.io>`_: Powerful and popular library for gradient boosted trees;
  includes native support for distributed training using dask.
- `lightgbm <https://lightgbm.readthedocs.io>`_: Similar to XGBoost; lightgmb also natively supplies native
  distributed training for decision trees.

Deploying Dask
--------------
There are many different implementations of the Dask distributed cluster.

- `dask-jobqueue <https://jobqueue.dask.org>`_: Deploy Dask on job queuing systems like PBS, Slurm, MOAB, SGE, LSF, and HTCondor.
- `dask-kubernetes <https://kubernetes.dask.org>`_: Deploy Dask workers on Kubernetes from within a Python script or interactive session.
- `dask-helm <https://helm.dask.org>`_: Deploy Dask and (optionally) Jupyter or JupyterHub on Kubernetes easily using Helm.
- `dask-yarn / Hadoop <https://yarn.dask.org>`_: Deploy Dask on YARN clusters, such as are found in traditional Hadoop
  installations.
- `dask-cloudprovider <https://cloudprovider.dask.org>`_: Deploy Dask on various cloud platforms such as AWS, Azure, and GCP
  leveraging cloud native APIs.
- `dask-gateway <https://gateway.dask.org>`_: Secure, multi-tenant server for managing Dask clusters. Launch and use Dask
  clusters in a shared, centrally managed cluster environment, without requiring users to have direct access to the underlying
  cluster backend.
- `dask-cuda <https://github.com/rapidsai/dask-cuda>`_: Construct a Dask cluster which resembles ``LocalCluster``  and is specifically
  optimized for GPUs.
Python API
==========

You can create a ``dask.distributed`` scheduler by importing and creating a
``Client`` with no arguments.  This overrides whatever default was previously
set.

.. code-block:: python

   from dask.distributed import Client
   client = Client()

You can navigate to ``http://localhost:8787/status`` to see the diagnostic
dashboard if you have Bokeh installed.

Client
------

You can trivially set up a local cluster on your machine by instantiating a Dask
Client with no arguments

.. code-block:: python

   from dask.distributed import Client
   client = Client()

This sets up a scheduler in your local process along with a number of workers and
threads per worker related to the number of cores in your machine.

If you want to run workers in your same process, you can pass the
``processes=False`` keyword argument.

.. code-block:: python

   client = Client(processes=False)

This is sometimes preferable if you want to avoid inter-worker communication
and your computations release the GIL.  This is common when primarily using
NumPy or Dask Array.


LocalCluster
------------

The ``Client()`` call described above is shorthand for creating a LocalCluster
and then passing that to your client.

.. code-block:: python

   from dask.distributed import Client, LocalCluster
   cluster = LocalCluster()
   client = Client(cluster)

This is equivalent, but somewhat more explicit.

You may want to look at the
keyword arguments available on ``LocalCluster`` to understand the options available
to you on handling the mixture of threads and processes, like specifying explicit
ports, and so on.

Cluster manager features
------------------------

Instantiating a cluster manager class like ``LocalCluster`` and then passing it to the
``Client`` is a common pattern. Cluster managers also provide useful utilities to help
you understand what is going on.

For example you can retreive the Dashboard URL.

.. code-block:: python

   >>> cluster.dashboard_link
   'http://127.0.0.1:8787/status'

You can retreive logs from cluster components.

.. code-block:: python

   >>> cluster.get_logs()
   {'Cluster': '',
   'Scheduler': "distributed.scheduler - INFO - Clear task state\ndistributed.scheduler - INFO -   S...

If you are using a cluster manager that supports scaling you can modify the number of workers manually
or automatically based on workload.

.. code-block:: python

   >>> cluster.scale(10)  # Sets the number of workers to 10

   >>> cluster.adapt(minimum=1, maximum=10)  # Allows the cluster to auto scale to 10 when tasks are computed

Reference
---------

.. currentmodule:: distributed.deploy.local

.. autoclass:: LocalCluster
   :members:
.. _order:

Ordering
========

.. note::

   This is an advanced topic that most users won't need to worry about.

When Dask is given a task graph to compute, it needs to choose an order to
execute the tasks in. We have some constraints: dependencies must be executed
before their dependants. But beyond that there's a large space of options. We
want Dask to choose an ordering that maximizes parallelism while minimizing
the footprint necessary to run a computation.

At a high level, Dask has a policy that works towards *small goals* with *big steps*.

1.  **Small goals**: prefer tasks that have few total dependents and whose final
    dependents have few total dependencies.

    We prefer to prioritize those tasks that help branches of computation that
    can terminate quickly.

    With more detail, we compute the total number of dependencies that each
    task depends on (both its own dependencies, the dependencies of its
    dependencies, and so on), and then we choose those tasks that drive towards
    results with a low number of total dependencies.  We choose to prioritize
    tasks that work towards finishing shorter computations first.

2.  **Big steps**: prefer tasks with many dependents

    However, many tasks work towards the same final dependents.  Among those,
    we choose those tasks with the most work left to do.  We want to finish
    the larger portions of a sub-computation before we start on the smaller
    ones.

This is done with :func:`dask.order.order`. A more technical discussion is available
in :ref:`scheduling-policy`. https://distributed.dask.org/en/latest/scheduling-policies.html
also discusses scheduling with a focus on the distributed scheduler, which includes
additional choices beyond the static ordering documented here.

Debugging
---------

Most of the time Dask's ordering does well. But this is a genuinely hard problem
and there might be cases where you observe unexpectedly high memory usage or
communication, which may be a result of poor ordering. This section describes
how you would identify an ordering problem, and some steps you can take to
mitigate the problem.

Consider a computation that loads several chains of data from disk independently,
stacks pieces of them together, and does some reduction:

.. code-block:: python

   >>> # create data on disk
   >>> import dask.array as da
   >>> x = da.zeros((12500, 10000), chunks=('10MB', -1))
   >>> da.to_zarr(x, 'saved_x1.zarr', overwrite=True)
   >>> da.to_zarr(x, 'saved_y1.zarr', overwrite=True)
   >>> da.to_zarr(x, 'saved_x2.zarr', overwrite=True)
   >>> da.to_zarr(x, 'saved_y2.zarr', overwrite=True)

We can load the data

.. code-block:: python

   >>> # load the data.
   >>> x1 = da.from_zarr('saved_x1.zarr')
   >>> y1 = da.from_zarr('saved_x2.zarr')
   >>> x2 = da.from_zarr('saved_y1.zarr')
   >>> y2 = da.from_zarr('saved_y2.zarr')

And do some computation on it

.. code-block:: python

   >>> def evaluate(x1, y1, x2, y2):
   ...     u = da.stack([x1, y1])
   ...     v = da.stack([x2, y2])
   ...     components = [u, v, u ** 2 + v ** 2]
   ...     return [
   ...         abs(c[0] - c[1]).mean(axis=-1)
   ...         for c in components
   ...     ]
   >>> results = evaluate(x1, y1, x2, y2)

You can use :func:`dask.visualize` with ``color="order"`` to visualize a
task graph with the static ordering included as node labels. As usual with
``dask.visualize``, you may need to trim down the problem to a smaller size,
so we'll slice off a subset of the data. Make sure to include ``optimize_graph=True``
to get a true representation of what order the tasks will be executed in.


.. code-block:: python

   >>> import dask
   >>> n = 125 * 4
   >>> dask.visualize(evaluate(x1[:n], y1[:n], x2[:n], y2[:n]),
   ...                optimize_graph=True, color="order",
   ...                cmap="autumn", node_attr={"penwidth": "4"})


.. image:: images/order-failure.png
   :alt: Complex task graph of several vertical node chains at the output, and a few input sub-trees. In between these sections, there is a many-to-many area of crossing dependency arrows. The color coding of the output trees is interleaved without a clear progression.

In this visualization the nodes are colored by order of execution (from dark red
to light yellow) and the node labels are the order Dask's assigned to each task.

It's a bit hard to see, but there are actually four mostly independent "towers"
of execution here. We start at the middle-right array (label 1, bottom), move
up to the right (label 8, top-right) and then jump to a completely different
array (label 11, bottom-left). However, computing the first tower (downstream
of label 8, top-right) required loading some data from our second input array
(label 5, bottom-right). We'd much prefer to finish tasks downstream of it.

When Dask is executing that task graph, you might observe high memory usage.
The poor static ordering means we fail to complete tasks that would let us
release pieces of data. We load more pieces into memory at once, leading to
higher memory usage.

This specific ordering failure (which may be fixed) comes from the shared
dependencies (the boxes at the bottom of each task, which represent the input
Zarr arrays) at the bottom of each computation chain. We can inline those and
see the effect of ordering:

.. code-block:: python

   >>> # load and profile data
   >>> x1 = da.from_zarr('saved_x1.zarr', inline_array=True)
   >>> y1 = da.from_zarr('saved_x2.zarr', inline_array=True)
   >>> x2 = da.from_zarr('saved_y1.zarr', inline_array=True)
   >>> y2 = da.from_zarr('saved_y2.zarr', inline_array=True)

   >>> import dask
   >>> n = 125 * 4
   >>> dask.visualize(evaluate(x1[:n], y1[:n], x2[:n], y2[:n]),
   ...                optimize_graph=True, color="order",
   ...                cmap="autumn", node_attr={"penwidth": "4"})


.. image:: images/order-success.png
   :alt: Complex task graph of several vertical node chains at the output, and a similar number of input blocks. The outputs and inputs are linked by simple nodes of a few inputs each, laid out without significant crossover between sections of the tree. The color coding of the output chains shows clear progression in the order of execution with each output color having a corresponding input of the same color.

At a glance, we can see that this ordering is looks much more regular and
uniform. There's fewer lines crossing, and the color of the ordering moves
smoothly from bottom to top, left to right. This shows that Dask is completing
one chain of computation before moving onto the next.

The lesson here is *not* "always use ``inline_array=True``". While the static
ordering looks better, there are other :ref:`phases-of-computation` to consider.
Whether the actual performance is better will depend on more factors than we
can consider here. See :func:`dask.array.from_array` for more.

Instead, the lessons to take away here are:

1. What symptoms might lead you to diagnose Dask's ordering as a problem (e.g.
   high memory usage)
2. How to generate and read task graphs with Dask's ordering information
   included.
Best Practices
==============

It is easy to get started with Dask arrays, but using them *well* does require
some experience.  This page contains suggestions for best practices, and
includes solutions to common problems.

Use NumPy
---------

If your data fits comfortably in RAM and you are not performance bound, then
using NumPy might be the right choice.  Dask adds another layer of complexity
which may get in the way.

If you are just looking for speedups rather than scalability then you may want
to consider a project like `Numba <https://numba.pydata.org>`_


Select a good chunk size
------------------------

A common performance problem among Dask Array users is that they have chosen a
chunk size that is either too small (leading to lots of overhead) or poorly
aligned with their data (leading to inefficient reading).

While optimal sizes and shapes are highly problem specific, it is rare to see
chunk sizes below 100 MB in size.  If you are dealing with float64 data then
this is around ``(4000, 4000)`` in size for a 2D array or ``(100, 400, 400)``
for a 3D array.

You want to choose a chunk size that is large in order to reduce the number of
chunks that Dask has to think about (which affects overhead) but also small
enough so that many of them can fit in memory at once.  Dask will often have as
many chunks in memory as twice the number of active threads.


Orient your chunks
------------------

When reading data you should align your chunks with your storage format.
Most array storage formats store data in chunks themselves.  If your Dask array
chunks aren't multiples of these chunk shapes then you will have to read the
same data repeatedly, which can be expensive.  Note though that often
storage formats choose chunk sizes that are much smaller than is ideal for
Dask, closer to 1MB than 100MB.  In these cases you should choose a Dask chunk
size that aligns with the storage chunk size and that every Dask chunk
dimension is a multiple of the storage chunk dimension.

So for example if we have an HDF file that has chunks of size ``(128, 64)``, we
might choose a chunk shape of ``(1280, 6400)``.

.. code-block:: python

   >>> import h5py
   >>> storage = h5py.File('myfile.hdf5')['x']
   >>> storage.chunks
   (128, 64)

   >>> import dask.array as da
   >>> x = da.from_array(storage, chunks=(1280, 6400))

Note that if you provide ``chunks='auto'`` then Dask Array will look for a
``.chunks`` attribute and use that to provide a good chunking.


Avoid Oversubscribing Threads
-----------------------------

By default Dask will run as many concurrent tasks as you have logical cores.
It assumes that each task will consume about one core.  However, many
array-computing libraries are themselves multi-threaded, which can cause
contention and low performance.  In particular the BLAS/LAPACK libraries that
back most of NumPy's linear algebra routines are often multi-threaded, and need
to be told to use only one thread explicitly.  You can do this with the
following environment variables (using bash ``export`` command below, but this
may vary depending on your operating system).

.. code-block:: bash

   export OMP_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   export OPENBLAS_NUM_THREADS=1

You need to run this before you start your Python process for it to take effect.


Consider Xarray
---------------

The `Xarray <http://xarray.pydata.org/en/stable/>`_ package wraps around Dask
Array, and so offers the same scalability, but also adds convenience when
dealing with complex datasets.  In particular Xarray can help with the
following:

1.  Manage multiple arrays together as a consistent dataset
2.  Read from a stack of HDF or NetCDF files at once
3.  Switch between Dask Array and NumPy with a consistent API

Xarray is used in wide range of fields, including physics, astronomy, geoscience,
microscopy, bioinformatics, engineering, finance, and deep learning.
Xarray also has a thriving user community that is good at providing support.


Build your own Operations
-------------------------

Often we want to perform computations for which there is no exact function
in Dask Array.  In these cases we may be able to use some of the more generic
functions to build our own.  These include:

.. currentmodule:: dask.array

.. autosummary::
   blockwise
   map_blocks
   map_overlap
   reduction

These functions may help you to apply a function that you write for NumPy
functions onto larger Dask arrays.
:orphan:

Diagnostics (distributed)
=========================

The :doc:`Dask distributed scheduler <scheduling>` provides live feedback in two
forms:

1.  An interactive dashboard containing many plots and tables with live
    information
2.  A progress bar suitable for interactive use in consoles or notebooks

Dashboard
---------

.. raw:: html

    <iframe width="560"
            height="315"
            src="https://www.youtube.com/embed/N_GqzcuGLCY"
            frameborder="0"
            allow="autoplay; encrypted-media"
            allowfullscreen>
    </iframe>

If `Bokeh <https://docs.bokeh.org>`_ is installed
then the dashboard will start up automatically whenever the scheduler is created.
For local use this happens when you create a client with no arguments:

.. code-block:: python

   from dask.distributed import Client
   client = Client()  # start distributed scheduler locally.  Launch dashboard

It is typically served at ``http://localhost:8787/status`` ,
but may be served elsewhere if this port is taken.
The address of the dashboard will be displayed if you are in a Jupyter Notebook,
or can be queried from ``client.dashboard_link``
(or for older versions of distributed, ``client.scheduler_info()['services']``).

There are numerous pages with information about task runtimes, communication,
statistical profiling, load balancing, memory use, and much more.
For more information we recommend the video guide above.

.. currentmodule:: dask.distributed

.. autosummary::
   Client


Capture diagnostics
-------------------

.. autosummary::
   get_task_stream
   Client.profile
   performance_report

You can capture some of the same information that the dashboard presents for
offline processing using the ``get_task_stream`` and ``Client.profile``
functions.  These capture the start and stop time of every task and transfer,
as well as the results of a statistical profiler.

.. code-block:: python

   with get_task_stream(plot='save', filename="task-stream.html") as ts:
       x.compute()

   client.profile(filename="dask-profile.html")

   history = ts.data

Additionally, Dask can save many diagnostics dashboards at once including the
task stream, worker profiles, bandwidths, etc. with the ``performance_report``
context manager:

.. code-block:: python

    from dask.distributed import performance_report

    with performance_report(filename="dask-report.html"):
        ## some dask computation

The following video demonstrates the ``performance_report`` context manager in greater
detail:

.. raw:: html

    <iframe width="560"
            height="315"
            src="https://www.youtube.com/embed/nTMGbkS761Q"
            frameborder="0"
            allow="autoplay; encrypted-media"
            allowfullscreen>
    </iframe>


Progress bar
------------

.. currentmodule:: dask.distributed

.. autosummary::
   progress

The ``dask.distributed`` progress bar differs from the ``ProgressBar`` used for
:doc:`local diagnostics <diagnostics-local>`.
The ``progress`` function takes a Dask object that is executing in the background:

.. code-block:: python

   # Progress bar on a single-machine scheduler
   from dask.diagnostics import ProgressBar

   with ProgressBar():
       x.compute()

   # Progress bar with the distributed scheduler
   from dask.distributed import Client, progress

   client = Client()  # use dask.distributed by default

   x = x.persist()  # start computation in the background
   progress(x)      # watch progress

   x.compute()      # convert to final result when done if desired


Connecting to the Dashboard
---------------------------

Some computer networks may restrict access to certain ports or only allow
access from certain machines.  If you are unable to access the dashboard then
you may want to contact your IT administrator.

Some common problems and solutions follow:

Specify an accessible port
~~~~~~~~~~~~~~~~~~~~~~~~~~

Some clusters restrict the ports that are visible to the outside world.  These
ports may include the default port for the web interface, ``8787``.  There are
a few ways to handle this:

1.  Open port ``8787`` to the outside world.  Often this involves asking your
    cluster administrator.
2.  Use a different port that is publicly accessible using the
    ``--dashboard-address :8787`` option on the ``dask-scheduler`` command.
3.  Use fancier techniques, like `Port Forwarding`_

Port Forwarding
~~~~~~~~~~~~~~~

If you have SSH access then one way to gain access to a blocked port is through
SSH port forwarding. A typical use case looks like the following:

.. code:: bash

   local$ ssh -L 8000:localhost:8787 user@remote
   remote$ dask-scheduler  # now, the web UI is visible at localhost:8000
   remote$ # continue to set up dask if needed -- add workers, etc

It is then possible to go to ``localhost:8000`` and see Dask Web UI. This same approach is
not specific to dask.distributed, but can be used by any service that operates over a
network, such as Jupyter notebooks. For example, if we chose to do this we could
forward port 8888 (the default Jupyter port) to port 8001 with
``ssh -L 8001:localhost:8888 user@remote``.

Required Packages
~~~~~~~~~~~~~~~~~

Bokeh must be installed in your scheduler's environment to run the dashboard. If it's not the dashboard page will instruct you to install it.

Depending on your configuration, you might also need to install ``jupyter-server-proxy`` to access the dashboard.

API
---

.. autofunction:: progress
.. autofunction:: get_task_stream
High Performance Computers
==========================

Relevant Machines
-----------------

This page includes instructions and guidelines when deploying Dask on high
performance supercomputers commonly found in scientific and industry research
labs.  These systems commonly have the following attributes:

1.  Some mechanism to launch MPI applications or use job schedulers like
    SLURM, SGE, TORQUE, LSF, DRMAA, PBS, or others
2.  A shared network file system visible to all machines in the cluster
3.  A high performance network interconnect, such as Infiniband
4.  Little or no node-local storage


Where to start
--------------

Most of this page documents various ways and best practices to use Dask on an
HPC cluster.  This is technical and aimed both at users with some experience
deploying Dask and also system administrators.

The preferred and simplest way to run Dask on HPC systems today both for new,
experienced users or administrator is to use
`dask-jobqueue <https://jobqueue.dask.org>`_.

However, dask-jobqueue is slightly oriented toward interactive analysis usage,
and it might be better to use tools like dask-mpi in some routine batch
production workloads.


Dask-jobqueue and Dask-drmaa
----------------------------

`dask-jobqueue <https://jobqueue.dask.org>`_ provides cluster managers for PBS,
SLURM, LSF, SGE and other resource managers. You can launch a Dask cluster on
these systems like this.

.. code-block:: python

   from dask_jobqueue import PBSCluster

   cluster = PBSCluster(cores=36,
                        memory="100GB",
                        project='P48500028',
                        queue='premium',
                        interface='ib0',
                        walltime='02:00:00')

   cluster.scale(100)  # Start 100 workers in 100 jobs that match the description above

   from dask.distributed import Client
   client = Client(cluster)    # Connect to that cluster

Dask-jobqueue provides a lot of possibilities like adaptive dynamic scaling
of workers, we recommend reading the `dask-jobqueue documentation
<https://jobqueue.dask.org>`_ first to get a basic system running and then
returning to this documentation for fine-tuning if necessary.


Using MPI
---------

You can launch a Dask cluster using ``mpirun`` or ``mpiexec`` and the
`dask-mpi <http://mpi.dask.org/en/latest/>`_ command line tool.

.. code-block:: bash

   mpirun --np 4 dask-mpi --scheduler-file /home/$USER/scheduler.json

.. code-block:: python

   from dask.distributed import Client
   client = Client(scheduler_file='/path/to/scheduler.json')

This depends on the `mpi4py <https://mpi4py.readthedocs.io/>`_ library.  It only
uses MPI to start the Dask cluster and not for inter-node communication. MPI
implementations differ: the use of ``mpirun --np 4`` is specific to the
``mpich`` or ``open-mpi`` MPI implementation installed through conda and linked
to mpi4py.

.. code-block:: bash

   conda install mpi4py

It is not necessary to use exactly this implementation, but you may want to
verify that your ``mpi4py`` Python library is linked against the proper
``mpirun/mpiexec`` executable and that the flags used (like ``--np 4``) are
correct for your system.  The system administrator of your cluster should be
very familiar with these concerns and able to help.

In some setups, MPI processes are not allowed to fork other processes. In this
case, we recommend using ``--no-nanny`` option in order to prevent dask from
using an additional nanny process to manage workers.

Run ``dask-mpi --help`` to see more options for the ``dask-mpi`` command.


Using a Shared Network File System and a Job Scheduler
------------------------------------------------------

.. note:: This section is not necessary if you use a tool like dask-jobqueue.

Some clusters benefit from a shared File System (NFS, GPFS, Lustre or alike),
and can use this to communicate the scheduler location to the workers::

   dask-scheduler --scheduler-file /path/to/scheduler.json  # writes address to file

   dask-worker --scheduler-file /path/to/scheduler.json  # reads file for address
   dask-worker --scheduler-file /path/to/scheduler.json  # reads file for address

.. code-block:: python

   >>> client = Client(scheduler_file='/path/to/scheduler.json')

This can be particularly useful when deploying ``dask-scheduler`` and
``dask-worker`` processes using a job scheduler like
SGE/SLURM/Torque/etc.  Here is an example using SGE's ``qsub`` command::

    # Start a dask-scheduler somewhere and write the connection information to a file
    qsub -b y /path/to/dask-scheduler --scheduler-file /home/$USER/scheduler.json

    # Start 100 dask-worker processes in an array job pointing to the same file
    qsub -b y -t 1-100 /path/to/dask-worker --scheduler-file /home/$USER/scheduler.json

Note, the ``--scheduler-file`` option is *only* valuable if your scheduler and
workers share a network file system.


High Performance Network
------------------------

Many HPC systems have both standard Ethernet networks as well as
high-performance networks capable of increased bandwidth.  You can instruct
Dask to use the high-performance network interface by using the ``--interface``
keyword with the ``dask-worker``, ``dask-scheduler``, or ``dask-mpi`` commands or
the ``interface=`` keyword with the dask-jobqueue ``Cluster`` objects:

.. code-block:: bash

   mpirun --np 4 dask-mpi --scheduler-file /home/$USER/scheduler.json --interface ib0

In the code example above, we have assumed that your cluster has an Infiniband
network interface called ``ib0``. You can check this by asking your system
administrator or by inspecting the output of ``ifconfig``

.. code-block:: bash

	$ ifconfig
	lo          Link encap:Local Loopback                       # Localhost
				inet addr:127.0.0.1  Mask:255.0.0.0
				inet6 addr: ::1/128 Scope:Host
	eth0        Link encap:Ethernet  HWaddr XX:XX:XX:XX:XX:XX   # Ethernet
				inet addr:192.168.0.101
				...
	ib0         Link encap:Infiniband                           # Fast InfiniBand
				inet addr:172.42.0.101

https://stackoverflow.com/questions/43881157/how-do-i-use-an-infiniband-network-with-dask


Local Storage
-------------

Users often exceed memory limits available to a specific Dask deployment.  In
normal operation, Dask spills excess data to disk, often to the default
temporary directory.

However, in HPC systems this default temporary directory may point to an
network file system (NFS) mount which can cause problems as Dask tries to read
and write many small files.  *Beware, reading and writing many tiny files from
many distributed processes is a good way to shut down a national
supercomputer*.

If available, it's good practice to point Dask workers to local storage, or
hard drives that are physically on each node.  Your IT administrators will be
able to point you to these locations.  You can do this with the
``--local-directory`` or ``local_directory=`` keyword in the ``dask-worker``
command::

   dask-mpi ... --local-directory /path/to/local/storage

or any of the other Dask Setup utilities, or by specifying the
following :doc:`configuration value <../../configuration>`:

.. code-block:: yaml

   temporary-directory: /path/to/local/storage

However, not all HPC systems have local storage.  If this is the case then you
may want to turn off Dask's ability to spill to disk altogether.  See `this
page <https://distributed.dask.org/en/latest/worker.html#memory-management>`_
for more information on Dask's memory policies.  Consider changing the
following values in your ``~/.config/dask/distributed.yaml`` file to disable
spilling data to disk:

.. code-block:: yaml

   distributed:
     worker:
       memory:
         target: false  # don't spill to disk
         spill: false  # don't spill to disk
         pause: 0.80  # pause execution at 80% memory use
         terminate: 0.95  # restart the worker at 95% use

This stops Dask workers from spilling to disk, and instead relies entirely on
mechanisms to stop them from processing when they reach memory limits.

As a reminder, you can set the memory limit for a worker using the
``--memory-limit`` keyword::

   dask-mpi ... --memory-limit 10GB


Launch Many Small Jobs
----------------------

.. note:: This section is not necessary if you use a tool like dask-jobqueue.

HPC job schedulers are optimized for large monolithic jobs with many nodes that
all need to run as a group at the same time.  Dask jobs can be quite a bit more
flexible: workers can come and go without strongly affecting the job.  If we
split our job into many smaller jobs, we can often get through the job
scheduling queue much more quickly than a typical job.  This is particularly
valuable when we want to get started right away and interact with a Jupyter
notebook session rather than waiting for hours for a suitable allocation block
to become free.

So, to get a large cluster quickly, we recommend allocating a dask-scheduler
process on one node with a modest wall time (the intended time of your session)
and then allocating many small single-node dask-worker jobs with shorter wall
times (perhaps 30 minutes) that can easily squeeze into extra space in the job
scheduler.  As you need more computation, you can add more of these single-node
jobs or let them expire.


Use Dask to co-launch a Jupyter server
--------------------------------------

Dask can help you by launching other services alongside it.  For example, you
can run a Jupyter notebook server on the machine running the ``dask-scheduler``
process with the following commands

.. code-block:: python

   from dask.distributed import Client
   client = Client(scheduler_file='scheduler.json')

   import socket
   host = client.run_on_scheduler(socket.gethostname)

   def start_jlab(dask_scheduler):
       import subprocess
       proc = subprocess.Popen(['/path/to/jupyter', 'lab', '--ip', host, '--no-browser'])
       dask_scheduler.jlab_proc = proc

   client.run_on_scheduler(start_jlab)
Kubernetes Native
=================

See external documentation on Dask-Kubernetes_ for more information.

.. _Dask-Kubernetes: https://kubernetes.dask.org
Kubernetes
==========

.. toctree::
   :maxdepth: 1
   :hidden:

   Helm <deploying-kubernetes-helm.rst>
   Native <deploying-kubernetes-native.rst>

Kubernetes_ is a popular system for deploying distributed applications on clusters,
particularly in the cloud.  You can use Kubernetes to launch Dask workers in the
following two ways:

1.  **Helm**:

    You can deploy Dask and (optionally) Jupyter or JupyterHub on Kubernetes
    easily using Helm_

    .. code-block:: bash

       helm repo add dask https://helm.dask.org/    # add the Dask Helm chart repository
       helm repo update                             # get latest Helm charts
       # For single-user deployments, use dask/dask
       helm install my-dask dask/dask               # deploy standard Dask chart
       # For multi-user deployments, use dask/daskhub
       helm install my-dask dask/daskhub            # deploy JupyterHub & Dask

    This is a good choice if you want to do the following:

    1.  Run a managed Dask cluster for a long period of time
    2.  Also deploy a Jupyter / JupyterHub server from which to run code
    3.  Share the same Dask cluster between many automated services
    4.  Try out Dask for the first time on a cloud-based system
        like Amazon, Google, or Microsoft Azure where you already have
        a Kubernetes cluster. If you don't already have Kubernetes deployed,
        see our :doc:`Cloud documentation <cloud>`.

    You can also use the ``HelmCluster`` cluster manager from dask-kubernetes to manage your
    Helm Dask cluster from within your Python session.

    .. code-block:: python

       from dask_kubernetes import HelmCluster

       cluster = HelmCluster(release_name="myrelease")
       cluster.scale(10)

    .. note::

      For more information, see :doc:`Dask and Helm documentation <kubernetes-helm>`.

2.  **Native**:
    You can quickly deploy Dask workers on Kubernetes
    from within a Python script or interactive session using Dask-Kubernetes_

    .. code-block:: python

       from dask_kubernetes import KubeCluster
       cluster = KubeCluster.from_yaml('worker-template.yaml')
       cluster.scale(20)  # add 20 workers
       cluster.adapt()    # or create and destroy workers dynamically based on workload

       from dask.distributed import Client
       client = Client(cluster)

    This is a good choice if you want to do the following:

    1.  Dynamically create a personal and ephemeral deployment for interactive use
    2.  Allow many individuals the ability to launch their own custom dask deployments,
        rather than depend on a centralized system
    3.  Quickly adapt Dask cluster size to the current workload

    .. note::

      For more information, see Dask-Kubernetes_ documentation.

You may also want to see the documentation on using
:doc:`Dask with Docker containers <docker>`
to help you manage your software environments on Kubernetes.

.. _Kubernetes: https://kubernetes.io/
.. _Dask-Kubernetes: https://kubernetes.dask.org/
.. _Helm: https://helm.sh/
Internal Design
===============

Overview
--------

.. image:: images/array.svg
   :width: 40 %
   :align: right
   :alt: 12 rectangular blocks arranged as a 4-row, 3-column layout. Each block includes 'x' and its location in the table starting with ('x',0,0) in the top-left, and a size of 5x8.

Dask arrays define a large array with a grid of blocks of smaller arrays.
These arrays may be actual arrays or functions that produce arrays. We 
define a Dask array with the following components:

*  A Dask graph with a special set of keys designating blocks
   such as ``('x', 0, 0), ('x', 0, 1), ...`` (See :doc:`Dask graph
   documentation <graphs>` for more details)
*  A sequence of chunk sizes along each dimension called ``chunks``,
   for example ``((5, 5, 5, 5), (8, 8, 8))``
*  A name to identify which keys in the Dask graph refer to this array, like
   ``'x'``
*  A NumPy dtype

Example
~~~~~~~

.. code-block:: python

   >>> import dask.array as da
   >>> x = da.arange(0, 15, chunks=(5,))

   >>> x.name
   'arange-539766a'

   >>> x.__dask_graph__()
   <dask.highlevelgraph.HighLevelGraph at 0x7f9f6f686d68>

   >>> dict(x.__dask_graph__())  # somewhat simplified
   {('arange-539766a', 0): (np.arange, 0, 5),
    ('arange-539766a', 1): (np.arange, 5, 10),
    ('arange-539766a', 2): (np.arange, 10, 15)}

   >>> x.chunks
   ((5, 5, 5),)

   >>> x.dtype
   dtype('int64')


Keys of the Dask graph
----------------------

By special convention, we refer to each block of the array with a tuple of the
form ``(name, i, j, k)``, with ``i, j, k`` being the indices of the block
ranging from ``0`` to the number of blocks in that dimension.  The Dask graph
must hold key-value pairs referring to these keys.  Moreover, it likely also
holds other key-value pairs required to eventually compute the desired values
(usually organised in a :doc:`HighLevelGraph <high-level-graphs>`, but shown
in a flattened form here for illustration):

.. code-block:: python

   {
    ('x', 0, 0): (add, 1, ('y', 0, 0)),
    ('x', 0, 1): (add, 1, ('y', 0, 1)),
    ...
    ('y', 0, 0): (getitem, dataset, (slice(0, 1000), slice(0, 1000))),
    ('y', 0, 1): (getitem, dataset, (slice(0, 1000), slice(1000, 2000)))
    ...
   }

The name of an ``Array`` object can be found in the ``name`` attribute.  One
can get a nested list of keys with the ``.__dask_keys__()`` method.  Additionally, 
one can flatten down this list with ``dask.array.core.flatten()``. This is sometimes
useful when building new dictionaries.

Chunks
------

We also store the size of each block along each axis.  This is composed of 
a tuple of tuples such that the length of the outer tuple is equal to the 
number of dimensions of the array, and the lengths of the inner tuples are 
equal to the number of blocks along each dimension.  In the example illustrated 
above this value is as follows::

    chunks = ((5, 5, 5, 5), (8, 8, 8))

Note that these numbers do not necessarily need to be regular.  We often create
regularly sized grids but blocks change shape after complex slicing.  Beware
that some operations do expect certain symmetries in the block-shapes.  For
example, matrix multiplication requires that blocks on each side have
anti-symmetric shapes.

Some ways in which ``chunks`` reflects properties of our array:

1.  ``len(x.chunks) == x.ndim``: the length of chunks is the number of dimensions
2.  ``tuple(map(sum, x.chunks)) == x.shape``: the sum of each internal chunk is the
    length of that dimension
3.  The length of each internal chunk is the number of keys in that dimension.
    For instance, for ``chunks == ((a, b), (d, e, f))`` and name == ``'x'``
    our array has tasks with the following keys::

       ('x', 0, 0), ('x', 0, 1), ('x', 0, 2)
       ('x', 1, 0), ('x', 1, 1), ('x', 1, 2)


Create an Array Object
----------------------

In order to create an ``da.Array`` object we need a graph with these special
keys::

    layer = {('x', 0, 0): ...}
    dsk = HighLevelGraph.from_collections('x', layer, dependencies=())

a name specifying which keys this array refers to::

    name = 'x'

and a chunks tuple::

    chunks = ((5, 5, 5, 5), (8, 8, 8))

Then, using these elements, one can construct an array::

    x = da.Array(dsk, name, chunks)

In short, ``dask.array`` operations update Dask graphs, update dtypes, and track chunk
shapes.


Example - ``eye`` function
--------------------------

As an example, let's build the ``np.eye`` function for ``dask.array`` to make the
identity matrix:

.. code-block:: python

   def eye(n, blocksize):
       chunks = ((blocksize,) * (n // blocksize),
                 (blocksize,) * (n // blocksize))

       name = 'eye' + next(tokens)  # unique identifier

       layer = {(name, i, j): (np.eye, blocksize)
                              if i == j else
                              (np.zeros, (blocksize, blocksize))
                for i in range(n // blocksize)
                for j in range(n // blocksize)}
       dsk = dask.highlevelgraph.HighLevelGraph.from_collections(name, layer, dependencies=())

       dtype = np.eye(0).dtype  # take dtype default from numpy

       return dask.array.Array(dsk, name, chunks, dtype)
Overlapping Computations
========================

Some array operations require communication of borders between neighboring
blocks.  Example operations include the following:

*  Convolve a filter across an image
*  Sliding sum/mean/max, ...
*  Search for image motifs like a Gaussian blob that might span the border of a
   block
*  Evaluate a partial derivative
*  Play the game of Life_

Dask Array supports these operations by creating a new array where each
block is slightly expanded by the borders of its neighbors.  This costs an
excess copy and the communication of many small chunks, but allows localized
functions to evaluate in an embarrassingly parallel manner.

The main API for these computations is the ``map_overlap`` method defined
below:

.. currentmodule:: dask.array

.. autosummary::
   map_overlap

.. autofunction:: map_overlap


Explanation
-----------

Consider two neighboring blocks in a Dask array:

.. image:: images/unoverlapping-neighbors.svg
   :width: 30%
   :alt: Two neighboring blocks which do not overlap.

We extend each block by trading thin nearby slices between arrays:

.. image:: images/overlapping-neighbors.svg
   :width: 30%
   :alt: Two neighboring block with thin strips along their shared border representing data shared between them.

We do this in all directions, including also diagonal interactions with the
overlap function:

.. image:: images/overlapping-blocks.svg
   :width: 40%
   :alt: A two-dimensional grid of blocks where each one has thin strips around their borders representing data shared from their neighbors. They include small corner bits for data shared from diagonal neighbors as well.

.. code-block:: python

   >>> import dask.array as da
   >>> import numpy as np

   >>> x = np.arange(64).reshape((8, 8))
   >>> d = da.from_array(x, chunks=(4, 4))
   >>> d.chunks
   ((4, 4), (4, 4))

   >>> g = da.overlap.overlap(d, depth={0: 2, 1: 1},
   ...                       boundary={0: 100, 1: 'reflect'})
   >>> g.chunks
   ((8, 8), (6, 6))

   >>> np.array(g)
   array([[100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
          [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
          [  0,   0,   1,   2,   3,   4,   3,   4,   5,   6,   7,   7],
          [  8,   8,   9,  10,  11,  12,  11,  12,  13,  14,  15,  15],
          [ 16,  16,  17,  18,  19,  20,  19,  20,  21,  22,  23,  23],
          [ 24,  24,  25,  26,  27,  28,  27,  28,  29,  30,  31,  31],
          [ 32,  32,  33,  34,  35,  36,  35,  36,  37,  38,  39,  39],
          [ 40,  40,  41,  42,  43,  44,  43,  44,  45,  46,  47,  47],
          [ 16,  16,  17,  18,  19,  20,  19,  20,  21,  22,  23,  23],
          [ 24,  24,  25,  26,  27,  28,  27,  28,  29,  30,  31,  31],
          [ 32,  32,  33,  34,  35,  36,  35,  36,  37,  38,  39,  39],
          [ 40,  40,  41,  42,  43,  44,  43,  44,  45,  46,  47,  47],
          [ 48,  48,  49,  50,  51,  52,  51,  52,  53,  54,  55,  55],
          [ 56,  56,  57,  58,  59,  60,  59,  60,  61,  62,  63,  63],
          [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100],
          [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]])


Boundaries
----------

With respect to overlapping, you can specify how to handle the boundaries.  Current policies
include the following:

*  ``periodic`` - wrap borders around to the other side
*  ``reflect`` - reflect each border outwards
*  ``any-constant`` - pad the border with this value

An example boundary kind argument might look like the following:

.. code-block:: python

   {0: 'periodic',
    1: 'reflect',
    2: np.nan}

Alternatively, you can use :py:func:`dask.array.pad` for other types of
paddings.


Map a function across blocks
----------------------------

Overlapping goes hand-in-hand with mapping a function across blocks.  This
function can now use the additional information copied over from the neighbors
that is not stored locally in each block:

.. code-block:: python

   >>> from scipy.ndimage.filters import gaussian_filter
   >>> def func(block):
   ...    return gaussian_filter(block, sigma=1)

   >>> filt = g.map_blocks(func)

While in this case we used a SciPy function, any arbitrary function could have been 
used instead. This is a good interaction point with Numba_.

If your function does not preserve the shape of the block, then you will need to
provide a ``chunks`` keyword argument. If your block size is regular, then this
argument can take a block shape of, for example, ``(1000, 1000)``. In case of
irregular block sizes, it must be a tuple with the full chunks shape like
``((1000, 700, 1000), (200, 300))``.

.. code-block:: python

   >>> g.map_blocks(myfunc, chunks=(5, 5))

If your function needs to know the location of the block on which it operates,
you can give your function a keyword argument ``block_id``:

.. code-block:: python

   def func(block, block_id=None):
       ...

This extra keyword argument will be given a tuple that provides the block
location like ``(0, 0)`` for the upper-left block or ``(0, 1)`` for the block
just to the right of that block.


Trim Excess
-----------

After mapping a blocked function, you may want to trim off the borders from each
block by the same amount by which they were expanded.  The function
``trim_internal`` is useful here and takes the same ``depth`` argument
given to ``overlap``:

.. code-block:: python

   >>> x.chunks
   ((10, 10, 10, 10), (10, 10, 10, 10))

   >>> y = da.overlap.trim_internal(x, {0: 2, 1: 1})
   >>> y.chunks
   ((6, 6, 6, 6), (8, 8, 8, 8))


Full Workflow
-------------

And so, a pretty typical overlapping workflow includes ``overlap``, ``map_blocks``
and ``trim_internal``:

.. code-block:: python

   >>> x = ...
   >>> g = da.overlap.overlap(x, depth={0: 2, 1: 2},
   ...                       boundary={0: 'periodic', 1: 'periodic'})
   >>> g2 = g.map_blocks(myfunc)
   >>> result = da.overlap.trim_internal(g2, {0: 2, 1: 2})

.. _Life: https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life
.. _Numba: https://numba.pydata.org/
GPUs
====

Dask works with GPUs in a few ways.


Custom Computations
-------------------

Many people use Dask alongside GPU-accelerated libraries like PyTorch and
TensorFlow to manage workloads across several machines.  They typically use
Dask's custom APIs, notably :doc:`Delayed <delayed>` and :doc:`Futures
<futures>`.

Dask doesn't need to know that these functions use GPUs.  It just runs Python
functions.  Whether or not those Python functions use a GPU is orthogonal to
Dask.  It will work regardless.

As a worked example, you may want to view this talk:

.. raw:: html

    <video width="560" height="315" controls>
        <source src="https://developer.download.nvidia.com/video/gputechconf/gtc/2019/video/S9198/s9198-dask-and-v100s-for-fast-distributed-batch-scoring-of-computer-vision-workloads.mp4"
                type="video/mp4">
    </video>

High Level Collections
----------------------

Dask can also help to scale out large array and dataframe computations by
combining the Dask Array and DataFrame collections with a GPU-accelerated
array or dataframe library.

Recall that :doc:`Dask Array <array>` creates a large array out of many NumPy
arrays and :doc:`Dask DataFrame <dataframe>` creates a large dataframe out of
many Pandas dataframes.  We can use these same systems with GPUs if we swap out
the NumPy/Pandas components with GPU-accelerated versions of those same
libraries, as long as the GPU accelerated version looks enough like
NumPy/Pandas in order to interoperate with Dask.

Fortunately, libraries that mimic NumPy, Pandas, and Scikit-Learn on the GPU do
exist.


DataFrames
~~~~~~~~~~

The `RAPIDS <https://rapids.ai>`_ libraries provide a GPU accelerated
Pandas-like library,
`cuDF <https://github.com/rapidsai/cudf>`_,
which interoperates well and is tested against Dask DataFrame.

If you have cuDF installed then you should be able to convert a Pandas-backed
Dask DataFrame to a cuDF-backed Dask DataFrame as follows:

.. code-block:: python

   import cudf

   df = df.map_partitions(cudf.from_pandas)  # convert pandas partitions into cudf partitions

However, cuDF does not support the entire Pandas interface, and so a variety of
Dask DataFrame operations will not function properly. Check the
`cuDF API Reference <https://docs.rapids.ai/api/cudf/stable/>`_
for currently supported interface.


Arrays
~~~~~~

.. note:: Dask's integration with CuPy relies on features recently added to
   NumPy and CuPy, particularly in version ``numpy>=1.17`` and ``cupy>=6``

`Chainer's CuPy <https://cupy.chainer.org/>`_ library provides a GPU
accelerated NumPy-like library that interoperates nicely with Dask Array.

If you have CuPy installed then you should be able to convert a NumPy-backed
Dask Array into a CuPy backed Dask Array as follows:

.. code-block:: python

   import cupy

   x = x.map_blocks(cupy.asarray)

CuPy is fairly mature and adheres closely to the NumPy API.  However, small
differences do exist and these can cause Dask Array operations to function
improperly. Check the
`CuPy Reference Manual <https://docs-cupy.chainer.org/en/stable/reference/index.html>`_
for API compatibility.


Scikit-Learn
~~~~~~~~~~~~

There are a variety of GPU accelerated machine learning libraries that follow
the Scikit-Learn Estimator API of fit, transform, and predict.  These can
generally be used within `Dask-ML's <https://ml.dask.org>`_ meta estimators,
such as `hyper parameter optimization <https://ml.dask.org/hyper-parameter-search.html>`_.

Some of these include:

-  `Skorch <https://skorch.readthedocs.io/>`_
-  `cuML <https://rapidsai.github.io/projects/cuml/en/latest/>`_
-  `LightGBM <https://github.com/Microsoft/LightGBM>`_
-  `XGBoost <https://xgboost.readthedocs.io/en/latest/>`_
-  `Thunder SVM <https://github.com/Xtra-Computing/thundersvm>`_
-  `Thunder GBM <https://github.com/Xtra-Computing/thundergbm>`_


Setup
-----

From the examples above we can see that the user experience of using Dask with
GPU-backed libraries isn't very different from using it with CPU-backed
libraries.  However, there are some changes you might consider making when
setting up your cluster.

Restricting Work
~~~~~~~~~~~~~~~~

By default Dask allows as many tasks as you have CPU cores to run concurrently.
However if your tasks primarily use a GPU then you probably want far fewer
tasks running at once.  There are a few ways to limit parallelism here:

-   Limit the number of threads explicitly on your workers using the
    ``--nthreads`` keyword in the CLI or the ``ncores=`` keyword the
    Cluster constructor.
-   Use `worker resources <https://distributed.dask.org/en/latest/resources.html>`_ and tag certain
    tasks as GPU tasks so that the scheduler will limit them, while leaving the
    rest of your CPU cores for other work

Specifying GPUs per Machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some configurations may have many GPU devices per node.  Dask is often used to
balance and coordinate work between these devices.

In these situations it is common to start one Dask worker per device, and use
the CUDA environment variable ``CUDA_VISIBLE_DEVICES`` to pin each worker to
prefer one device.

.. code-block:: bash

   # If we have four GPUs on one machine
   CUDA_VISIBLE_DEVICES=0 dask-worker ...
   CUDA_VISIBLE_DEVICES=1 dask-worker ...
   CUDA_VISIBLE_DEVICES=2 dask-worker ...
   CUDA_VISIBLE_DEVICES=3 dask-worker ...

The `Dask CUDA <https://github.com/rapidsai/dask-cuda>`_ project contains some
convenience CLI and Python utilities to automate this process.

Work in Progress
----------------

GPU computing is a quickly moving field today and as a result the information
in this page is likely to go out of date quickly.  We encourage interested
readers to check out `Dask's Blog <https://blog.dask.org>`_ which has more
timely updates on ongoing work.
.. _array.assignment:

Assignment
==========

Dask Array supports most of the NumPy assignment indexing syntax. In
particular, it supports combinations of the following:

* Indexing by integers: ``x[1] = y``
* Indexing by slices: ``x[2::-1] = y``
* Indexing by a list of integers: ``x[[0, -1, 1]] = y``
* Indexing by a 1-d :class:`numpy` array of integers: ``x[np.arange(3)] = y``
* Indexing by a 1-d :class:`~dask.array.Array` of integers: ``x[da.arange(3)] = y``, ``x[da.from_array([0, -1, 1])] = y``, ``x[da.where(np.array([1, 2, 3]) < 3)[0]] = y``
* Indexing by a list of booleans: ``x[[False, True, True]] = y``
* Indexing by a 1-d :class:`numpy` array of booleans: ``x[np.arange(3) > 0] = y``

It also supports:

* Indexing by one broadcastable :class:`~dask.array.Array` of
  booleans: ``x[x > 0] = y``.

However, it does not currently support the following:

* Indexing with lists in multiple axes: ``x[[1, 2, 3], [3, 1, 2]] = y``


.. _array.assignment.broadcasting:

Broadcasting
------------

The normal NumPy broadcasting rules apply:

.. code-block:: python

   >>> x = da.zeros((2, 6))
   >>> x[0] = 1
   >>> x[..., 1] = 2.0
   >>> x[:, 2] = [3, 4]
   >>> x[:, 5:2:-2] = [[6, 5]]
   >>> x.compute()
   array([[1., 2., 3., 5., 1., 6.],
          [0., 2., 4., 5., 0., 6.]])
   >>> x[1] = -x[0]
   >>> x.compute()
   array([[ 1.,  2.,  3.,  5.,  1.,  6.],
          [-1., -2., -3., -5., -1., -6.]])

.. _array.assignment.masking:

Masking
-------

Elements may be masked by assigning to the NumPy masked value, or to an
array with masked values:

.. code-block:: python

   >>> x = da.ones((2, 6))
   >>> x[0, [1, -2]] = np.ma.masked
   >>> x[1] = np.ma.array([0, 1, 2, 3, 4, 5], mask=[0, 1, 1, 0, 0, 0])
   >>> print(x.compute())
   [[1.0 -- 1.0 1.0 -- 1.0]
    [0.0 -- -- 3.0 4.0 5.0]]
   >>> x[:, 0] = x[:, 1]
   >>> print(x.compute())
   [[1.0 -- 1.0 1.0 -- 1.0]
    [0.0 -- -- 3.0 4.0 5.0]]
   >>> x[:, 0] = x[:, 1]
   >>> print(x.compute())
   [[-- -- 1.0 1.0 -- 1.0]
    [-- -- -- 3.0 4.0 5.0]]

If, and only if, a single broadcastable :class:`~dask.array.Array` of
booleans is provided then masked array assignment does not yet work as
expected. In this case the data underlying the mask are assigned:

.. code-block:: python

   >>> x = da.arange(12).reshape(2, 6)
   >>> x[x > 7] = np.ma.array(-99, mask=True)
   >>> print(x.compute())
   [[  0   1   2   3   4   5]
    [  6   7 -99 -99 -99 -99]]

Note that masked assignments do work when a boolean
:class:`~dask.array.Array` index used in a tuple, or implicit tuple,
of indices:

.. code-block:: python

   >>> x = da.arange(12).reshape(2, 6)
   >>> x[1, x[0] > 3] = np.ma.masked
   >>> print(x.compute())
   [[0 1 2 3 4 5]
    [6 7 8 9 -- --]]
   >>> x = da.arange(12).reshape(2, 6)
   >>> print(x.compute())
   [[ 0  1  2  3  4  5]
    [ 6  7  8  9 10 11]]
   >>> x[(x[:, 2] < 4,)] = np.ma.masked
   >>> print(x.compute())
   [[-- -- -- -- -- --]
    [6 7 8 9 10 11]]


Joins
=====

DataFrame joins are a common and expensive computation that benefit from a
variety of optimizations in different situations.  Understanding how your data
is laid out and what you're trying to accomplish can have a large impact on
performance.  This documentation page goes through the various different
options and their performance impacts.

Large to Large Unsorted Joins
-----------------------------

In the worst case scenario you have two large tables with many partitions each
and you want to join them both along a column that may not be sorted.

This can be slow.  In this case Dask DataFrame will need to move all of your
data around so that rows with matching values in the joining columns are in the
same partition.  This large-scale movement can create communication costs, and
can require a large amount of memory.  If enough memory can not be found then
Dask will have to read and write data to disk, which may cause other
performance costs.

These problems are solvable, but will be significantly slower than many other
operations.  They are best avoided if possible.

Large to Small Joins
--------------------

Many join or merge computations combine a large table with one small one.  If
the small table is either a single partition Dask DataFrame or even just a
normal Pandas DataFrame then the computation can proceed in an embarrassingly
parallel way, where each partition of the large DataFrame is joined against the
single small table.  This incurs almost no overhead relative to Pandas joins.

If your smaller table can easily fit in memory, then you might want to ensure
that it is a single partition with the repartition method.

.. code-block:: python

    import dask
    large = dask.datasets.timeseries(freq="10s", npartitions=10)
    small = dask.datasets.timeseries(freq="1D", dtypes={"z": int})

    small = small.repartition(npartitions=1)
    result = large.merge(small, how="left", on=["timestamp"])

Sorted Joins
------------

The Pandas merge API supports the ``left_index=`` and ``right_index=`` options
to perform joins on the index.  For Dask DataFrames these keyword options hold
special significance if the index has known divisions
(see :ref:`dataframe-design-partitions`).
In this case the DataFrame partitions are aligned along these divisions (which
is generally fast) and then an embarrassingly parallel Pandas join happens
across partition pairs.  This is generally relatively fast.

Sorted or indexed joins are a good solution to the large-large join problem.
If you plan to join against a dataset repeatedly then it may be worthwhile to
set the index ahead of time, and possibly store the data in a format that
maintains that index, like Parquet.

.. code-block:: python

    import dask
    import dask.dataframe as dd

    left = dask.datasets.timeseries(dtypes={"foo": int})

    # timeseries returns a dataframe indexed by
    # timestamp, we don't need to set_index.

    # left.set_index("timestamp")

    left.to_parquet("left", overwrite=True)
    left = dd.read_parquet("left")

    # If the dataframe can fit in RAM, you can also use persist

    # left = left.persist()

    right_one = dask.datasets.timeseries(dtypes={"bar": int})
    right_two = dask.datasets.timeseries(dtypes={"baz": int})

    result = left.merge(
        right_one, how="left", left_index=True, right_index=True)
    result = result.merge(
        right_two, how="left", left_index=True, right_index=True)
.. _graphs:

Task Graphs
===========

Internally, Dask encodes algorithms in a simple format involving Python dicts,
tuples, and functions. This graph format can be used in isolation from the
dask collections. Working directly with dask graphs is rare, though, unless you intend
to develop new modules with Dask.  Even then, :doc:`dask.delayed <delayed>` is
often a better choice. If you are a *core developer*, then you should start here.

.. toctree::
   :maxdepth: 1

   spec.rst
   custom-graphs.rst
   optimize.rst
   graph_manipulation.rst
   custom-collections.rst
   high-level-graphs.rst


Motivation
----------

Normally, humans write programs and then compilers/interpreters interpret them
(for example, ``python``, ``javac``, ``clang``).  Sometimes humans disagree with how
these compilers/interpreters choose to interpret and execute their programs.
In these cases, humans often bring the analysis, optimization, and execution of
code into the code itself.

Commonly a desire for parallel execution causes this shift of responsibility
from compiler to human developer.  In these cases, we often represent the
structure of our program explicitly as data within the program itself.

A common approach to parallel execution in user-space is *task scheduling*.  In
task scheduling we break our program into many medium-sized tasks or units of
computation, often a function call on a non-trivial amount of data.  We
represent these tasks as nodes in a graph with edges between nodes if one task
depends on data produced by another.  We call upon a *task scheduler* to
execute this graph in a way that respects these data dependencies and leverages
parallelism where possible, so multiple independent tasks can be run
simultaneously.

|

.. figure:: images/map-reduce-task-scheduling.svg
   :scale: 40%

   There are a number of methods for task scheduling, including embarrassingly parallel, MapReduce, and full task scheduling.

|

Many solutions exist.  This is a common approach in parallel execution
frameworks.  Often task scheduling logic hides within other larger frameworks
(e.g. Luigi, Storm, Spark, IPython Parallel, etc.) and so is often reinvented.
Dask is a specification that encodes full task scheduling with minimal incidental
complexity using terms common to all Python projects, namely, dicts, tuples,
and callables.  Ideally this minimum solution is easy to adopt and understand
by a broad community.

Example
-------

Consider the following simple program:

.. code-block:: python

   def inc(i):
       return i + 1

   def add(a, b):
       return a + b

   x = 1
   y = inc(x)
   z = add(y, 10)

We encode this as a dictionary in the following way:

.. code-block:: python

   d = {'x': 1,
        'y': (inc, 'x'),
        'z': (add, 'y', 10)}

Which is represented by the following Dask graph:

.. image:: _static/dask-simple.png
   :height: 400px
   :alt: A simple dask dictionary

|

While less pleasant than our original code, this representation can be analyzed
and executed by other Python code, not just the CPython interpreter.  We don't
recommend that users write code in this way, but rather that it is an
appropriate target for automated systems.  Also, in non-toy examples, the
execution times are likely much larger than for ``inc`` and ``add``, warranting
the extra complexity.


Schedulers
----------

The Dask library currently contains a few schedulers to execute these
graphs.  Each scheduler works differently, providing different performance
guarantees and operating in different contexts.  These implementations are not
special and others can write different schedulers better suited to other
applications or architectures easily.  Systems that emit dask graphs (like
Dask Array, Dask Bag, and so on) may leverage the appropriate scheduler for
the application and hardware.


Task Expectations
-----------------

When a task is submitted to Dask for execution, there are a number of assumptions
that are made about that task.

Don't Modify Data In-Place
~~~~~~~~~~~~~~~~~~~~~~~~~~

In general, tasks with side-effects that alter the state of a future in-place
are not recommended. Modifying data that is stored in Dask in-place can have
unintended consequences. For example, consider a workflow involving a Numpy
array:

.. code-block:: python

   from dask.distributed import Client
   import numpy as np

   client = Client()
   x = client.submit(np.arange, 10)  # [0, 1, 2, 3, ...]

   def f(arr):
       arr[arr > 5] = 0  # modifies input directly without making a copy
       arr += 1          # modifies input directly without making a copy
       return arr

   y = client.submit(f, x)

In the example above Dask will update the values of the Numpy array
``x`` in-place.  While efficient, this behavior can have unintended consequences,
particularly if other tasks need to use ``x``, or if Dask needs to rerun this
computation multiple times because of worker failure.


Avoid Holding the GIL
~~~~~~~~~~~~~~~~~~~~~

Some Python functions that wrap external C/C++ code can hold onto the GIL,
which stops other Python code from running in the background.  This is
troublesome because while Dask workers run your function, they also need to
communicate to each other in the background.

If you wrap external code then please try to release the GIL.  This is usually
easy to do if you are using any of the common solutions to code-wrapping like
Cython, Numba, ctypes or others.
Opportunistic Caching
=====================

Dask usually removes intermediate values as quickly as possible in order to
make space for more data to flow through your computation.  However, in some
cases, we may want to hold onto intermediate values, because they might be
useful for future computations in an interactive session.

We need to balance the following concerns:

1.  Intermediate results might be useful in future unknown computations
2.  Intermediate results also fill up memory, reducing space for the rest of our
    current computation

Negotiating between these two concerns helps us to leverage the memory that we
have available to speed up future, unanticipated computations.  Which intermediate results
should we keep?

This document explains an experimental, opportunistic caching mechanism that automatically
picks out and stores useful tasks.


Motivating Example
------------------

Consider computing the maximum value of a column in a CSV file:

.. code-block:: python

   >>> import dask.dataframe as dd
   >>> df = dd.read_csv('myfile.csv')
   >>> df.columns
   ['first-name', 'last-name', 'amount', 'id', 'timestamp']

   >>> df.amount.max().compute()
   1000

Even though our full dataset may be too large to fit in memory, the single
``df.amount`` column may be small enough to hold in memory just in case it
might be useful in the future.  This is often the case during data exploration,
because we investigate the same subset of our data repeatedly before moving on.

For example, we may now want to find the minimum of the amount column:

.. code-block:: python

   >>> df.amount.min().compute()
   -1000

Under normal operations, this would need to read through the entire CSV file over
again.  This is somewhat wasteful and stymies interactive data exploration.


Two Simple Solutions
--------------------

If we know ahead of time that we want both the maximum and minimum, we can
compute them simultaneously.  Dask will share intermediates intelligently,
reading through the dataset only once:

.. code-block:: python

   >>> dd.compute(df.amount.max(), df.amount.min())
   (1000, -1000)

If we know that this column fits in memory, then we can also explicitly
compute the column and then continue forward with straight Pandas:

.. code-block:: python

   >>> amount = df.amount.compute()
   >>> amount.max()
   1000
   >>> amount.min()
   -1000

If either of these solutions work for you, great.  Otherwise, continue on for a third approach.


Automatic Opportunistic Caching
-------------------------------

Another approach is to watch *all* intermediate computations, and *guess* which
ones might be valuable to keep for the future.  Dask has an *opportunistic
caching mechanism* that stores intermediate tasks that show the following
characteristics:

1.  Expensive to compute
2.  Cheap to store
3.  Frequently used

We can activate a fixed sized cache as a callback_:

.. _callback: diagnostics-local.html#custom-callbacks

.. code-block:: python

   >>> from dask.cache import Cache
   >>> cache = Cache(2e9)  # Leverage two gigabytes of memory
   >>> cache.register()    # Turn cache on globally

Now the cache will watch every small part of the computation and judge the
value of that part based on the three characteristics listed above (expensive
to compute, cheap to store, and frequently used).

Dask will hold on to 2GB of the
best intermediate results it can find, evicting older results as better results
come in.  If the ``df.amount`` column fits in 2GB, then probably all of it will
be stored while we keep working on it.

If we start work on something else,
then the ``df.amount`` column will likely be evicted to make space for other
more timely results:

.. code-block:: python

   >>> df.amount.max().compute()  # slow the first time
   1000
   >>> df.amount.min().compute()  # fast because df.amount is in the cache
   -1000
   >>> df.id.nunique().compute()  # starts to push out df.amount from cache


Cache tasks, not expressions
----------------------------

This caching happens at the low-level scheduling layer, not the high-level
Dask DataFrame or Dask Array layer.  We don't explicitly cache the column
``df.amount``.  Instead, we cache the hundreds of small pieces of that column
that form the dask graph.  It could be that we end up caching only a fraction
of the column.

This means that the opportunistic caching mechanism described above works for *all* Dask
computations, as long as those computations employ a consistent naming scheme
(as all of Dask DataFrame, Dask Array, and Dask Delayed do).

You can see which tasks are held by the cache by inspecting the following
attributes of the cache object:

.. code-block:: python

   >>> cache.cache.data
   <stored values>
   >>> cache.cache.heap.heap
   <scores of items in cache>
   >>> cache.cache.nbytes
   <number of bytes per item in cache>

The cache object is powered by cachey_, a tiny library for opportunistic
caching.

.. _cachey: https://github.com/blaze/cachey

Disclaimer
----------
Opportunistic caching is not available when using the distributed scheduler.

Restricting your cache to a fixed size like 2GB requires Dask to accurately count
the size of each of our objects in memory.  This can be tricky, particularly
for Pythonic objects like lists and tuples, and for DataFrames that contain
object dtypes.

It is entirely possible that the caching mechanism will
*undercount* the size of objects, causing it to use up more memory than
anticipated, which can lead to blowing up RAM and crashing your session.

Development Guidelines
======================

Dask is a community maintained project.  We welcome contributions in the form
of bug reports, documentation, code, design proposals, and more.
This page provides resources on how best to contribute.

.. note:: Dask strives to be a welcoming community of individuals with diverse
   backgrounds. For more information on our values, please see our
   `code of conduct
   <https://github.com/dask/governance/blob/main/code-of-conduct.md>`_
   and
   `diversity statement <https://github.com/dask/governance/blob/main/diversity.md>`_

Where to ask for help
---------------------

Dask conversation happens in the following places:

#.  `Dask Discourse forum`_: for usage questions and general discussion
#.  `Stack Overflow #dask tag`_: for usage questions
#.  `GitHub Issue Tracker`_: for discussions around new features or established bugs
#.  `Dask Community Slack`_: for real-time discussion

For usage questions and bug reports we prefer the use of Discourse, Stack Overflow
and GitHub issues over Slack chat.  Discourse, GitHub and Stack Overflow are more easily
searchable by future users, so conversations had there can be useful to many more people
than just those directly involved.

.. _`Dask Discourse forum`: https://dask.discourse.group
.. _`Stack Overflow  #dask tag`: https://stackoverflow.com/questions/tagged/dask
.. _`GitHub Issue Tracker`: https://github.com/dask/dask/issues/
.. _`Dask Community Slack`: https://join.slack.com/t/dask/shared_invite/zt-mfmh7quc-nIrXL6ocgiUH2haLYA914g


Separate Code Repositories
--------------------------

Dask maintains code and documentation in a few git repositories hosted on the
GitHub ``dask`` organization, https://github.com/dask.  This includes the primary
repository and several other repositories for different components.  A
non-exhaustive list follows:

*  https://github.com/dask/dask: The main code repository holding parallel
   algorithms, the single-machine scheduler, and most documentation
*  https://github.com/dask/distributed: The distributed memory scheduler
*  https://github.com/dask/dask-ml: Machine learning algorithms
*  https://github.com/dask/s3fs: S3 Filesystem interface
*  https://github.com/dask/gcsfs: GCS Filesystem interface
*  https://github.com/dask/hdfs3: Hadoop Filesystem interface
*  ...

Git and GitHub can be challenging at first.  Fortunately good materials exist
on the internet.  Rather than repeat these materials here, we refer you to
Pandas' documentation and links on this subject at
https://pandas.pydata.org/pandas-docs/stable/contributing.html


Issues
------

The community discusses and tracks known bugs and potential features in the
`GitHub Issue Tracker`_.  If you have a new idea or have identified a bug, then
you should raise it there to start public discussion.

If you are looking for an introductory issue to get started with development,
then check out the `"good first issue" label`_, which contains issues that are good
for starting developers.  Generally, familiarity with Python, NumPy, Pandas, and
some parallel computing are assumed.

.. _`"good first issue" label`: https://github.com/dask/dask/labels/good%20first%20issue


Development Environment
-----------------------

Download code
~~~~~~~~~~~~~

Make a fork of the main `Dask repository <https://github.com/dask/dask>`_ and
clone the fork::

   git clone https://github.com/<your-github-username>/dask.git
   cd dask

You should also pull the latest git tags (this ensures ``pip``'s dependency resolver
can successfully install Dask)::

   git remote add upstream https://github.com/dask/dask.git
   git pull upstream main --tags

Contributions to Dask can then be made by submitting pull requests on GitHub.


Install
~~~~~~~

From the top level of your cloned Dask repository you can install a
local version of Dask, along with all necessary dependencies, using
pip or conda_

.. _conda: https://conda.io/

``pip``::

  python -m pip install -e ".[complete,test]"

``conda``::

  conda env create -n dask-dev -f continuous_integration/environment-3.9.yaml
  conda activate dask-dev
  python -m pip install --no-deps -e .


Run Tests
~~~~~~~~~

Dask uses py.test_ for testing.  You can run tests from the main dask directory
as follows::

   py.test dask --verbose --doctest-modules

.. _py.test: https://docs.pytest.org/en/latest/


Contributing to Code
--------------------

Dask maintains development standards that are similar to most PyData projects.  These standards include
language support, testing, documentation, and style.

Python Versions
~~~~~~~~~~~~~~~

Dask supports Python versions 3.7, 3.8, and 3.9.
Name changes are handled by the :file:`dask/compatibility.py` file.

Test
~~~~

Dask employs extensive unit tests to ensure correctness of code both for today
and for the future.  Test coverage is expected for all code contributions.

Tests are written in a py.test style with bare functions:

.. code-block:: python

   def test_fibonacci():
       assert fib(0) == 0
       assert fib(1) == 0
       assert fib(10) == 55
       assert fib(8) == fib(7) + fib(6)

       for x in [-3, 'cat', 1.5]:
           with pytest.raises(ValueError):
               fib(x)

These tests should compromise well between covering all branches and fail cases
and running quickly (slow test suites get run less often).

You can run tests locally by running ``py.test`` in the local dask directory::

   py.test dask

You can also test certain modules or individual tests for faster response::

   py.test dask/dataframe

   py.test dask/dataframe/tests/test_dataframe.py::test_rename_index

If you want the tests to run faster, you can run them in parallel using
``pytest-xdist``::

   py.test dask -n auto

Tests run automatically on the Travis.ci and Appveyor continuous testing
frameworks on every push to every pull request on GitHub.

Tests are organized within the various modules' subdirectories::

    dask/array/tests/test_*.py
    dask/bag/tests/test_*.py
    dask/bytes/tests/test_*.py
    dask/dataframe/tests/test_*.py
    dask/diagnostics/tests/test_*.py

For the Dask collections like Dask Array and Dask DataFrame, behavior is
typically tested directly against the NumPy or Pandas libraries using the
``assert_eq`` functions:

.. code-block:: python

   import numpy as np
   import dask.array as da
   from dask.array.utils import assert_eq

   def test_aggregations():
       nx = np.random.random(100)
       dx = da.from_array(nx, chunks=(10,))

       assert_eq(nx.sum(), dx.sum())
       assert_eq(nx.min(), dx.min())
       assert_eq(nx.max(), dx.max())
       ...

This technique helps to ensure compatibility with upstream libraries and tends
to be simpler than testing correctness directly.  Additionally, by passing Dask
collections directly to the ``assert_eq`` function rather than call compute
manually, the testing suite is able to run a number of checks on the lazy
collections themselves.


Docstrings
~~~~~~~~~~

User facing functions should roughly follow the numpydoc_ standard, including
sections for ``Parameters``, ``Examples``, and general explanatory prose.

By default, examples will be doc-tested.  Reproducible examples in documentation
is valuable both for testing and, more importantly, for communication of common
usage to the user.  Documentation trumps testing in this case and clear
examples should take precedence over using the docstring as testing space.
To skip a test in the examples add the comment ``# doctest: +SKIP`` directly
after the line.

.. code-block:: python

   def fib(i):
       """ A single line with a brief explanation

       A more thorough description of the function, consisting of multiple
       lines or paragraphs.

       Parameters
       ----------
       i: int
            A short description of the argument if not immediately clear

       Examples
       --------
       >>> fib(4)
       3
       >>> fib(5)
       5
       >>> fib(6)
       8
       >>> fib(-1)  # Robust to bad inputs
       ValueError(...)
       """

.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard

Docstrings are tested under Python 3.8 on GitHub Actions. You can test
docstrings with pytest as follows::

   py.test dask --doctest-modules

Docstring testing requires ``graphviz`` to be installed. This can be done via::

   conda install -y graphviz


Code Formatting
~~~~~~~~~~~~~~~

Dask uses several code linters (flake8, black, isort, pyupgrade), which are enforced by
CI. Developers should run them locally before they submit a PR, through the single
command ``pre-commit run --all-files``. This makes sure that linter versions and options
are aligned for all developers.

Optionally, you may wish to setup the `pre-commit hooks <https://pre-commit.com/>`_ to
run automatically when you make a git commit. This can be done by running::

   pre-commit install

from the root of the Dask repository. Now the code linters will be run each time you
commit changes. You can skip these checks with ``git commit --no-verify`` or with the
short version ``git commit -n``.


Contributing to Documentation
-----------------------------

Dask uses Sphinx_ for documentation, hosted on https://readthedocs.org .
Documentation is maintained in the RestructuredText markup language (``.rst``
files) in ``dask/docs/source``.  The documentation consists both of prose
and API documentation.

The documentation is automatically built, and a live preview is available,
for each pull request submitted to Dask. Additionally, you may also
build the documentation yourself locally by following the instructions outlined
below.

How to build the Dask documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To build the documentation locally, make a fork of the main
`Dask repository <https://github.com/dask/dask>`_, clone the fork::

  git clone https://github.com/<your-github-username>/dask.git
  cd dask/docs

Install the packages in ``requirements-docs.txt``.

Optionally create and activate a ``conda`` environment first::

  conda create -n daskdocs -c conda-forge python=3.8
  conda activate daskdocs

Install the dependencies with ``pip``::

  python -m pip install -r requirements-docs.txt

Then build the documentation with ``make``::

   make html

The resulting HTML files end up in the ``build/html`` directory.

You can now make edits to rst files and run ``make html`` again to update
the affected pages.


Dask CI Infrastructure
----------------------

Github Actions
~~~~~~~~~~~~~~

Dask uses Github Actions for Continuous Integration (CI) testing for each PR.
These CI builds will run the test suite across a variety of Python versions, operating
systems, and package dependency versions.  Addtionally, if a commit message
includes the phrase ``test-upstream``, then an additional CI build will be
triggered which uses the development versions of several dependencies
including: NumPy, pandas, fsspec, etc.

The CI workflows for Github Actions are defined in
`.github/workflows <https://github.com/dask/dask/tree/main/.github/workflows>`_
with additonal scripts and metadata located in `continuous_integration
<https://github.com/dask/dask/tree/main/continuous_integration>`_


GPU CI
~~~~~~

Pull requests are also tested with a GPU enabled CI environment provided by
NVIDIA: `gpuCI <https://gpuci.gpuopenanalytics.com/>`_.
Unlike Github Actions, the CI environment for gpuCI is controlled with the
`rapidsai/dask-build-environment <https://github.com/rapidsai/dask-build-environment/>`_
docker image.  When making commits to the
`dask-build-environment repo <https://github.com/rapidsai/dask-build-environment/>`_ , a new image is built.
The docker image building process can be monitored
`here <https://gpuci.gpuopenanalytics.com/job/dask/job/dask-build-environment/job/branch/job/dask-build-env-main/>`_.
Note, the ``dask-build-environment`` has two separate Dockerfiles for Dask
and Distributed similiarlly, gpuCI will run for both `Dask
<https://gpuci.gpuopenanalytics.com/job/dask/job/dask/job/prb/job/dask-prb/>`_
and `Distributed
<https://gpuci.gpuopenanalytics.com/job/dask/job/distributed/job/prb/job/distributed-prb/>`_

For each PR, gpuCI will run all tests decorated with the pytest marker
``@pytest.mark.gpu``.  This is configured in the `gpuci folder
<https://github.com/dask/dask/tree/main/continuous_integration/gpuci>`_ .
Like Github Actions, gpuCI will not run when first time contributors to Dask or
Distributed submit PRs.  In this case, the gpuCI bot will comment on the PR:

.. note:: Can one of the admins verify this patch?

.. image:: images/gputester-msg.png
   :alt: "Screenshot of a GitHub comment left by the GPUtester bot, where the comment says 'Can one of the admins verify this patch?'."

Dask Maintainers can then approve gpuCI builds for these PRs with following choices:

- To only approve the PR contributor for the current PR, leave a comment which states ``ok to test``
- To approve the current PR and all future PRs from the contributor, leave a comment which states ``add to allowlist``

For more information about gpuCI please consult the `docs page
<https://docs.rapids.ai/gpuci>`_


.. _Sphinx: https://www.sphinx-doc.org/
API
===

The ``dask.delayed`` interface consists of one function, ``delayed``:

- ``delayed`` wraps functions

   Wraps functions. Can be used as a decorator, or around function calls
   directly (i.e. ``delayed(foo)(a, b, c)``). Outputs from functions wrapped in
   ``delayed`` are proxy objects of type ``Delayed`` that contain a graph of
   all operations done to get to this result.

- ``delayed`` wraps objects

   Wraps objects. Used to create ``Delayed`` proxies directly.

``Delayed`` objects can be thought of as representing a key in the dask task
graph. A ``Delayed`` supports *most* python operations, each of which creates
another ``Delayed`` representing the result:

- Most operators (``*``, ``-``, and so on)
- Item access and slicing (``a[0]``)
- Attribute access (``a.size``)
- Method calls (``a.index(0)``)

Operations that aren't supported include:

- Mutating operators (``a += 1``)
- Mutating magics such as ``__setitem__``/``__setattr__`` (``a[0] = 1``, ``a.foo = 1``)
- Iteration. (``for i in a: ...``)
- Use as a predicate (``if a: ...``)

The last two points in particular mean that ``Delayed`` objects cannot be used for
control flow, meaning that no ``Delayed`` can appear in a loop or if statement.
In other words you can't iterate over a ``Delayed`` object, or use it as part of
a condition in an if statement, but ``Delayed`` object can be used in a body of a loop
or if statement (i.e. the example above is fine, but if ``data`` was a ``Delayed``
object it wouldn't be).
Even with this limitation, many workflows can easily be parallelized.

.. currentmodule:: dask.delayed

.. autosummary::
   delayed
   Delayed

.. autofunction:: delayed
.. autoclass:: Delayed
Configuration
=============

Taking full advantage of Dask sometimes requires user configuration.
This might be to control logging verbosity, specify cluster configuration,
provide credentials for security, or any of several other options that arise in
production.

Configuration is specified in one of the following ways:

1.  YAML files in ``~/.config/dask/`` or ``/etc/dask/``
2.  Environment variables like ``DASK_DISTRIBUTED__SCHEDULER__WORK_STEALING=True``
3.  Default settings within sub-libraries

This combination makes it easy to specify configuration in a variety of
settings ranging from personal workstations, to IT-mandated configuration, to
docker images.


Access Configuration
--------------------

.. currentmodule:: dask

.. autosummary::
   dask.config.get

Dask's configuration system is usually accessed using the ``dask.config.get`` function.
You can use ``.`` for nested access, for example:

.. code-block:: python

   >>> import dask
   >>> import dask.distributed  # populate config with distributed defaults

   >>> dask.config.get("distributed.client") # use `.` for nested access
   {'heartbeat': '5s', 'scheduler-info-interval': '2s'}

   >>> dask.config.get("distributed.scheduler.unknown-task-duration")
   '500ms'

You may wish to inspect the ``dask.config.config`` dictionary to get a sense
for what configuration is being used by your current system.

Note that the ``get`` function treats underscores and hyphens identically.
For example, ``dask.config.get("temporary-directory")`` is equivalent to
``dask.config.get("temporary_directory")``.

Values like ``"128 MiB"`` and ``"10s"`` are parsed using the functions in
:ref:`api.utilities`.

Specify Configuration
---------------------

YAML files
~~~~~~~~~~

You can specify configuration values in YAML files. For example:

.. code-block:: yaml

   array:
     chunk-size: 128 MiB

   distributed:
     worker:
       memory:
         spill: 0.85  # default: 0.7
         target: 0.75  # default: 0.6
         terminate: 0.98  # default: 0.95
            
     dashboard:
       # Locate the dashboard if working on a Jupyter Hub server
       link: /user/<user>/proxy/8787/status
        

These files can live in any of the following locations:

1.  The ``~/.config/dask`` directory in the user's home directory
2.  The ``{sys.prefix}/etc/dask`` directory local to Python
3.  The ``{prefix}/etc/dask`` directories with ``{prefix}`` in `site.PREFIXES
    <https://docs.python.org/3/library/site.html#site.PREFIXES>`_
4.  The root directory (specified by the ``DASK_ROOT_CONFIG`` environment
    variable or ``/etc/dask/`` by default)

Dask searches for *all* YAML files within each of these directories and merges
them together, preferring configuration files closer to the user over system
configuration files (preference follows the order in the list above).
Additionally, users can specify a path with the ``DASK_CONFIG`` environment
variable, which takes precedence at the top of the list above.

The contents of these YAML files are merged together, allowing different
Dask subprojects like ``dask-kubernetes`` or ``dask-ml`` to manage configuration
files separately, but have them merge into the same global configuration.


Environment Variables
~~~~~~~~~~~~~~~~~~~~~

You can also specify configuration values with environment variables like
the following:

.. code-block:: bash

   export DASK_DISTRIBUTED__SCHEDULER__WORK_STEALING=True
   export DASK_DISTRIBUTED__SCHEDULER__ALLOWED_FAILURES=5
   export DASK_DISTRIBUTED__DASHBOARD__LINK="/user/<user>/proxy/8787/status"

resulting in configuration values like the following:

.. code-block:: python

   {
       'distributed': {
           'scheduler': {
               'work-stealing': True,
               'allowed-failures': 5
           }
       }
   }

Dask searches for all environment variables that start with ``DASK_``, then
transforms keys by converting to lower case and changing double-underscores to
nested structures.

Dask tries to parse all values with `ast.literal_eval
<https://docs.python.org/3/library/ast.html#ast.literal_eval>`_, letting users
pass numeric and boolean values (such as ``True`` in the example above) as well
as lists, dictionaries, and so on with normal Python syntax.

Environment variables take precedence over configuration values found in YAML
files.

Defaults
~~~~~~~~

Additionally, individual subprojects may add their own default values when they
are imported.  These are always added with lower priority than the YAML files
or environment variables mentioned above:

.. code-block:: python

   >>> import dask.config
   >>> dask.config.config  # no configuration by default
   {}

   >>> import dask.distributed
   >>> dask.config.config  # New values have been added
   {
       'scheduler': ...,
       'worker': ...,
       'tls': ...
   }


Directly within Python
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   dask.config.set

Configuration is stored within a normal Python dictionary in
``dask.config.config`` and can be modified using normal Python operations.

Additionally, you can temporarily set a configuration value using the
``dask.config.set`` function.  This function accepts a dictionary as an input
and interprets ``"."`` as nested access:

.. code-block:: python

   >>> dask.config.set({'scheduler.work-stealing': True})

This function can also be used as a context manager for consistent cleanup:

.. code-block:: python

   with dask.config.set({'scheduler.work-stealing': True}):
       ...

Note that the ``set`` function treats underscores and hyphens identically.
For example, ``dask.config.set({'scheduler.work-stealing': True})`` is
equivalent to ``dask.config.set({'scheduler.work_stealing': True})``.

Distributing configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~

It may also be desirable to package up your whole Dask configuration for use on
another machine. This is used in some Dask Distributed libraries to ensure remote components
have the same configuration as your local system.

This is typically handled by the downstream libraries which use base64 encoding to pass
config via the ``DASK_INTERNAL_INHERIT_CONFIG`` environment variable.

.. autosummary::
   dask.config.serialize
   dask.config.deserialize

Conversion Utility
~~~~~~~~~~~~~~~~~~

It is possible to configure Dask inline with dot notation, with YAML or via environment variables. You can enter
your own configuration items below to convert back and forth.

.. warning::
   This utility is designed to improve understanding of converting between different notations
   and does not claim to be a perfect implementation. Please use for reference only.

**YAML**

.. raw:: html

   <textarea id="configConvertUtilYAML" name="configConvertUtilYAML" rows="10" cols="50" class="configTextArea" wrap="off">
   array:
      chunk-size: 128 MiB

   distributed:
      workers:
         memory:
            spill: 0.85
            target: 0.75
            terminate: 0.98
   </textarea>

**Environment variable**

.. raw:: html

   <textarea id="configConvertUtilEnv" name="configConvertUtilEnv" rows="10" cols="50" class="configTextArea" wrap="off">
   export DASK_ARRAY__CHUNK_SIZE="128 MiB"
   export DASK_DISTRIBUTED__WORKERS__MEMORY__SPILL=0.85
   export DASK_DISTRIBUTED__WORKERS__MEMORY__TARGET=0.75
   export DASK_DISTRIBUTED__WORKERS__MEMORY__TERMINATE=0.98
   </textarea>

**Inline with dot notation**

.. raw:: html

   <textarea id="configConvertUtilCode" name="configConvertUtilCode" rows="10" cols="50" class="configTextArea" wrap="off">
   >>> dask.config.set({"array.chunk-size": "128 MiB"})
   >>> dask.config.set({"distributed.workers.memory.spill": 0.85})
   >>> dask.config.set({"distributed.workers.memory.target": 0.75})
   >>> dask.config.set({"distributed.workers.memory.terminate": 0.98})
   </textarea>

Updating Configuration
----------------------

Manipulating configuration dictionaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   dask.config.merge
   dask.config.update
   dask.config.expand_environment_variables

As described above, configuration can come from many places, including several
YAML files, environment variables, and project defaults.  Each of these
provides a configuration that is possibly nested like the following:

.. code-block:: python

   x = {'a': 0, 'c': {'d': 4}}
   y = {'a': 1, 'b': 2, 'c': {'e': 5}}

Dask will merge these configurations respecting nested data structures, and
respecting order:

.. code-block:: python

   >>> dask.config.merge(x, y)
   {'a': 1, 'b': 2, 'c': {'d': 4, 'e': 5}}

You can also use the ``update`` function to update the existing configuration
in place with a new configuration.  This can be done with priority being given
to either config.  This is often used to update the global configuration in
``dask.config.config``:

.. code-block:: python

   dask.config.update(dask.config, new, priority='new')  # Give priority to new values
   dask.config.update(dask.config, new, priority='old')  # Give priority to old values

Sometimes it is useful to expand environment variables stored within a
configuration.  This can be done with the ``expand_environment_variables``
function:

.. code-block:: python

    dask.config.config = dask.config.expand_environment_variables(dask.config.config)

Refreshing Configuration
~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   dask.config.collect
   dask.config.refresh

If you change your environment variables or YAML files, Dask will not
immediately see the changes.  Instead, you can call ``refresh`` to go through
the configuration collection process and update the default configuration:

.. code-block:: python

   >>> dask.config.config
   {}

   >>> # make some changes to yaml files

   >>> dask.config.refresh()
   >>> dask.config.config
   {...}

This function uses ``dask.config.collect``, which returns the configuration
without modifying the global configuration.  You might use this to determine
the configuration of particular paths not yet on the config path:

.. code-block:: python

   >>> dask.config.collect(paths=[...])
   {...}

Downstream Libraries
--------------------

.. autosummary::
   dask.config.ensure_file
   dask.config.update
   dask.config.update_defaults

Downstream Dask libraries often follow a standard convention to use the central
Dask configuration.  This section provides recommendations for integration
using a fictional project, ``dask-foo``, as an example.

Downstream projects typically follow the following convention:

1.  Maintain default configuration in a YAML file within their source
    directory::

       setup.py
       dask_foo/__init__.py
       dask_foo/config.py
       dask_foo/core.py
       dask_foo/foo.yaml  # <---

2.  Place configuration in that file within a namespace for the project:

    .. code-block:: yaml

       # dask_foo/foo.yaml

       foo:
         color: red
         admin:
           a: 1
           b: 2

3.  Within a config.py file (or anywhere) load that default config file and
    update it into the global configuration:

    .. code-block:: python

       # dask_foo/config.py
       import os
       import yaml

       import dask.config

       fn = os.path.join(os.path.dirname(__file__), 'foo.yaml')

       with open(fn) as f:
           defaults = yaml.safe_load(f)

       dask.config.update_defaults(defaults)

4.  Within that same config.py file, copy the ``'foo.yaml'`` file to the user's
    configuration directory if it doesn't already exist.

    We also comment the file to make it easier for us to change defaults in the
    future.

    .. code-block:: python

       # ... continued from above

       dask.config.ensure_file(source=fn, comment=True)

    The user can investigate ``~/.config/dask/*.yaml`` to see all of the
    commented out configuration files to which they have access.

5.  Ensure that this file is run on import by including it in ``__init__.py``:

    .. code-block:: python

       # dask_foo/__init__.py

       from . import config

6.  Within ``dask_foo`` code, use the ``dask.config.get`` function to access
    configuration values:

    .. code-block:: python

       # dask_foo/core.py

       def process(fn, color=dask.config.get('foo.color')):
           ...

7.  You may also want to ensure that your yaml configuration files are included
    in your package.  This can be accomplished by including the following line
    in your MANIFEST.in::

       recursive-include <PACKAGE_NAME> *.yaml

    and the following in your setup.py ``setup`` call:

    .. code-block:: python

        from setuptools import setup

        setup(...,
              include_package_data=True,
              ...)

This process keeps configuration in a central place, but also keeps it safe
within namespaces.  It places config files in an easy to access location
by default (``~/.config/dask/\*.yaml``), so that users can easily discover what
they can change, but maintains the actual defaults within the source code, so
that they more closely track changes in the library.

However, downstream libraries may choose alternative solutions, such as
isolating their configuration within their library, rather than using the
global dask.config system.  All functions in the ``dask.config`` module also
work with parameters, and do not need to mutate global state.


API
---

.. autofunction:: dask.config.get
.. autofunction:: dask.config.set
.. autofunction:: dask.config.merge
.. autofunction:: dask.config.update
.. autofunction:: dask.config.collect
.. autofunction:: dask.config.refresh
.. autofunction:: dask.config.ensure_file
.. autofunction:: dask.config.expand_environment_variables


Configuration Reference
-----------------------

.. contents:: :local:

.. note::
   It is possible to configure Dask inline with dot notation, with YAML or via environment variables.
   See the `conversion utility <#conversion-utility>`_ for converting the following dot notation to other forms.

Dask
~~~~

.. dask-config-block::
    :location: dask
    :config: https://raw.githubusercontent.com/dask/dask/main/dask/dask.yaml
    :schema: https://raw.githubusercontent.com/dask/dask/main/dask/dask-schema.yaml


Distributed Client
~~~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.client
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml

Distributed Comm
~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.comm
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml


Distributed Dashboard
~~~~~~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.dashboard
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml


Distributed Deploy
~~~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.deploy
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml


Distributed Scheduler
~~~~~~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.scheduler
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml


Distributed Worker
~~~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.worker
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml


Distributed Nanny
~~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.nanny
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml


Distributed Admin
~~~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.admin
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml


Distributed RMM
~~~~~~~~~~~~~~~

.. dask-config-block::
    :location: distributed.rmm
    :config: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed.yaml
    :schema: https://raw.githubusercontent.com/dask/distributed/main/distributed/distributed-schema.yaml
Slicing
=======

Dask Array supports most of the NumPy slicing syntax.  In particular, it
supports the following:

*  Slicing by integers and slices: ``x[0, :5]``
*  Slicing by lists/arrays of integers: ``x[[1, 2, 4]]``
*  Slicing by lists/arrays of booleans: ``x[[False, True, True, False, True]]``
*  Slicing one :class:`~dask.array.Array` with an :class:`~dask.array.Array` of bools: ``x[x > 0]``
*  Slicing one :class:`~dask.array.Array` with a zero or one-dimensional :class:`~dask.array.Array`
   of ints: ``a[b.argtopk(5)]``

However, it does not currently support the following:

*  Slicing with lists in multiple axes: ``x[[1, 2, 3], [3, 2, 1]]``

   This is straightforward to add though.  If you have a use case then raise an
   issue. Also, users interested in this should take a look at
   :attr:`~dask.array.Array.vindex`.

*  Slicing one :class:`~dask.array.Array` with a multi-dimensional :class:`~dask.array.Array` of ints

.. _array.slicing.efficiency:

Efficiency
----------

The normal Dask schedulers are smart enough to compute only those blocks that
are necessary to achieve the desired slicing.  Hence, large operations may be cheap
if only a small output is desired.

In the example below, we create a Dask array with a trillion elements with million
element sized blocks.  We then operate on the entire array and finally slice out
only a portion of the output:

.. code-block:: python

   >>> # Trillion element array of ones, in 1000 by 1000 blocks
   >>> x = da.ones((1000000, 1000000), chunks=(1000, 1000))

   >>> da.exp(x)[:1500, :1500]
   ...

This only needs to compute the top-left four blocks to achieve the result.  We
are slightly wasteful on those blocks where we need only partial results.  Moreover,
we are also a bit wasteful in that we still need to manipulate the Dask graph
with a million or so tasks in it.  This can cause an interactive overhead of a
second or two.

Slicing with concrete indexers (a list of integers, say) has a couple of possible
failure modes that are worth mentioning. First, when you're indexing a chunked
axis, Dask will typically "match" the chunking on the output.

.. code-block:: python

   # Array of ones, chunked along axis 0
   >>> a = da.ones((4, 10000, 10000), chunks=(1, -1, -1))

If we slice that with a *sorted* sequence of integers, Dask will return one chunk
per input chunk (notice the output `chunksize` is 1, since the indices ``0``
and ``1`` are in separate chunks in the input).

   >>> a[[0, 1], :, :]          #doctest: +SKIP
   dask.array<getitem, shape=(2, 10000, 10000), dtype=float64, chunksize=(1, 10000, 10000), chunktype=numpy.ndarray>

But what about repeated indices? Dask continues to return one chunk per input chunk,
but if you have many repetitions from the same input chunk, your output chunk could
be much larger.

.. code-block:: python

   >>> a[[0] * 15, :, :]
   PerformanceWarning: Slicing is producing a large chunk. To accept the large
   chunk and silence this warning, set the option
       >>> with dask.config.set({'array.slicing.split_large_chunks': False}):
       ...     array[indexer]

   To avoid creating the large chunks, set the option
       >>> with dask.config.set({'array.slicing.split_large_chunks': True}):
       ...     array[indexer]
   dask.array<getitem, shape=(15, 10000, 10000), dtype=float64, chunksize=(15, 10000, 10000), chunktype=numpy.ndarray>

Previously we had a chunksize of ``1`` along the first dimension since we selected
just one element from each input chunk. But now we've selected 15 elements
from the first chunk, producing a large output chunk.

Dask warns when indexing like this produces a chunk that's 5x larger
than the ``array.chunk-size`` config option. You have two options to deal with
that warning:

1. Set ``dask.config.set({"array.slicing.split_large_chunks": False})`` to
   allow the large chunk and silence the warning.
2. Set ``dask.config.set({"array.slicing.split_large_chunks": True})`` to
   avoid creating the large chunk in the first place.

The right choice will depend on your downstream operations. See :ref:`array.chunks`
for more on choosing chunk sizes.
Scheduling
==========

All of the large-scale Dask collections like
:doc:`Dask Array <array>`, :doc:`Dask DataFrame <dataframe>`, and :doc:`Dask Bag <bag>`
and the fine-grained APIs like :doc:`delayed <delayed>` and :doc:`futures <futures>`
generate task graphs where each node in the graph is a normal Python function
and edges between nodes are normal Python objects
that are created by one task as outputs and used as inputs in another task.
After Dask generates these task graphs, it needs to execute them on parallel hardware.
This is the job of a *task scheduler*.
Different task schedulers exist, and each will consume a task graph and compute the
same result, but with different performance characteristics.

Dask has two families of task schedulers:

1.  **Single-machine scheduler**: This scheduler provides basic features on a
    local process or thread pool.  This scheduler was made first and is the
    default.  It is simple and cheap to use, although it can only be used on
    a single machine and does not scale
2.  **Distributed scheduler**: This scheduler is more sophisticated, offers
    more features, but also requires a bit more effort to set up.  It can
    run locally or distributed across a cluster

|

.. image:: images/dask-overview-schedulers.svg
   :alt: Dask is composed of three parts. "Collections" create "Task Graphs" which are then sent to the "Scheduler" for execution. There are two types of schedulers that are described in more detail below.
   :align: center
   :scale: 135%

|

For different computations you may find better performance with particular scheduler settings.
This document helps you understand how to choose between and configure different schedulers,
and provides guidelines on when one might be more appropriate.


Local Threads
-------------

.. code-block:: python

   import dask
   dask.config.set(scheduler='threads')  # overwrite default with threaded scheduler

The threaded scheduler executes computations with a local
``concurrent.futures.ThreadPoolExecutor``.
It is lightweight and requires no setup.
It introduces very little task overhead (around 50us per task)
and, because everything occurs in the same process,
it incurs no costs to transfer data between tasks.
However, due to Python's Global Interpreter Lock (GIL),
this scheduler only provides parallelism when your computation is dominated by non-Python code,
as is primarily the case when operating on numeric data in NumPy arrays, Pandas DataFrames,
or using any of the other C/C++/Cython based projects in the ecosystem.

The threaded scheduler is the default choice for
:doc:`Dask Array <array>`, :doc:`Dask DataFrame <dataframe>`, and :doc:`Dask Delayed <delayed>`.
However, if your computation is dominated by processing pure Python objects
like strings, dicts, or lists,
then you may want to try one of the process-based schedulers below
(we currently recommend the distributed scheduler on a local machine).


Local Processes
---------------

.. note::

   The distributed scheduler described a couple sections below is often a better choice today.
   We encourage readers to continue reading after this section.

.. code-block:: python

   import dask
   dask.config.set(scheduler='processes')  # overwrite default with multiprocessing scheduler


The multiprocessing scheduler executes computations with a local
``concurrent.futures.ProcessPoolExecutor``.
It is lightweight to use and requires no setup.
Every task and all of its dependencies are shipped to a local process,
executed, and then their result is shipped back to the main process.
This means that it is able to bypass issues with the GIL and provide parallelism
even on computations that are dominated by pure Python code,
such as those that process strings, dicts, and lists.

However, moving data to remote processes and back can introduce performance penalties,
particularly when the data being transferred between processes is large.
The multiprocessing scheduler is an excellent choice when workflows are relatively linear,
and so does not involve significant inter-task data transfer
as well as when inputs and outputs are both small, like filenames and counts.

This is common in basic data ingestion workloads,
such as those are common in :doc:`Dask Bag <bag>`,
where the multiprocessing scheduler is the default:

.. code-block:: python

   >>> import dask.bag as db
   >>> db.read_text('*.json').map(json.loads).pluck('name').frequencies().compute()
   {'alice': 100, 'bob': 200, 'charlie': 300}

For more complex workloads,
where large intermediate results may be depended upon by multiple downstream tasks,
we generally recommend the use of the distributed scheduler on a local machine.
The distributed scheduler is more intelligent about moving around large intermediate results.

.. _single-threaded-scheduler:

Single Thread
-------------

.. code-block:: python

   import dask
   dask.config.set(scheduler='synchronous')  # overwrite default with single-threaded scheduler

The single-threaded synchronous scheduler executes all computations in the local thread
with no parallelism at all.
This is particularly valuable for debugging and profiling,
which are more difficult when using threads or processes.

For example, when using IPython or Jupyter notebooks, the ``%debug``, ``%pdb``, or ``%prun`` magics
will not work well when using the parallel Dask schedulers
(they were not designed to be used in a parallel computing context).
However, if you run into an exception and want to step into the debugger,
you may wish to rerun your computation under the single-threaded scheduler
where these tools will function properly.


Dask Distributed (local)
------------------------

.. code-block:: python

   from dask.distributed import Client
   client = Client()
   # or
   client = Client(processes=False)

The Dask distributed scheduler can either be :doc:`setup on a cluster <how-to/deploy-dask-clusters>`
or run locally on a personal machine.  Despite having the name "distributed",
it is often pragmatic on local machines for a few reasons:

1.  It provides access to asynchronous API, notably :doc:`Futures <futures>`
2.  It provides a diagnostic dashboard that can provide valuable insight on
    performance and progress
3.  It handles data locality with more sophistication, and so can be more
    efficient than the multiprocessing scheduler on workloads that require
    multiple processes

You can read more about using the Dask distributed scheduler on a single machine in
:doc:`these docs <how-to/deploy-dask/single-distributed>`.


Dask Distributed (Cluster)
--------------------------

You can also run Dask on a distributed cluster.
There are a variety of ways to set this up depending on your cluster.
We recommend referring to :doc:`how to deploy Dask clusters <how-to/deploy-dask-clusters>` for more information.

.. _scheduling-configuration:

Configuration
-------------

You can configure the global default scheduler by using the ``dask.config.set(scheduler...)`` command.
This can be done globally:

.. code-block:: python

   dask.config.set(scheduler='threads')

   x.compute()

or as a context manager:

.. code-block:: python

   with dask.config.set(scheduler='threads'):
       x.compute()

or within a single compute call:

.. code-block:: python

   x.compute(scheduler='threads')

Each scheduler may support extra keywords specific to that scheduler. For example,
the pool-based single-machine scheduler allows you to provide custom pools or
specify the desired number of workers:

.. code-block:: python

   from concurrent.futures import ThreadPoolExecutor
   with dask.config.set(pool=ThreadPoolExecutor(4)):
       x.compute()

   with dask.config.set(num_workers=4):
       x.compute()

Note that Dask also supports custom ``concurrent.futures.Executor`` subclasses,
such as the ``ReusablePoolExecutor`` from loky_:

.. _loky: https://github.com/joblib/loky

.. code-block:: python

   from loky import get_reusable_executor
   with dask.config.set(scheduler=get_reusable_executor()):
       x.compute()

Other libraries like ipyparallel_ and mpi4py_ also supply
``concurrent.futures.Executor`` subclasses that could be used as well.

.. _ipyparallel: https://ipyparallel.readthedocs.io/en/latest/examples/Futures.html#Executors
.. _mpi4py: https://mpi4py.readthedocs.io/en/latest/mpi4py.futures.html
API
---

.. currentmodule:: dask.array

Top level functions
~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   add
   all
   allclose
   angle
   any
   append
   apply_along_axis
   apply_over_axes
   arange
   arccos
   arccosh
   arcsin
   arcsinh
   arctan
   arctan2
   arctanh
   argmax
   argmin
   argtopk
   argwhere
   around
   array
   asanyarray
   asarray
   atleast_1d
   atleast_2d
   atleast_3d
   average
   bincount
   bitwise_and
   bitwise_not
   bitwise_or
   bitwise_xor
   block
   blockwise
   broadcast_arrays
   broadcast_to
   cbrt
   coarsen
   ceil
   choose
   clip
   compress
   concatenate
   conj
   copysign
   corrcoef
   cos
   cosh
   count_nonzero
   cov
   cumprod
   cumsum
   deg2rad
   degrees
   diag
   diagonal
   diff
   divmod
   digitize
   dot
   dstack
   ediff1d
   einsum
   empty
   empty_like
   equal
   exp
   exp2
   expm1
   eye
   fabs
   fix
   flatnonzero
   flip
   flipud
   fliplr
   float_power
   floor
   floor_divide
   fmax
   fmin
   fmod
   frexp
   fromfunction
   frompyfunc
   full
   full_like
   gradient
   greater
   greater_equal
   histogram
   histogram2d
   histogramdd
   hstack
   hypot
   imag
   indices
   insert
   invert
   isclose
   iscomplex
   isfinite
   isin
   isinf
   isneginf
   isnan
   isnull
   isposinf
   isreal
   ldexp
   less
   linspace
   log
   log10
   log1p
   log2
   logaddexp
   logaddexp2
   logical_and
   logical_not
   logical_or
   logical_xor
   map_overlap
   map_blocks
   matmul
   max
   maximum
   mean
   median
   meshgrid
   min
   minimum
   mod
   modf
   moment
   moveaxis
   multiply
   nanargmax
   nanargmin
   nancumprod
   nancumsum
   nanmax
   nanmean
   nanmedian
   nanmin
   nanprod
   nanstd
   nansum
   nanvar
   nan_to_num
   negative
   nextafter
   nonzero
   not_equal
   notnull
   ones
   ones_like
   outer
   pad
   percentile
   ~core.PerformanceWarning
   piecewise
   power
   prod
   ptp
   rad2deg
   radians
   ravel
   real
   reciprocal
   rechunk
   reduction
   register_chunk_type
   remainder
   repeat
   reshape
   result_type
   rint
   roll
   rollaxis
   rot90
   round
   searchsorted
   sign
   signbit
   sin
   sinc
   sinh
   sqrt
   square
   squeeze
   stack
   std
   subtract
   sum
   take
   tan
   tanh
   tensordot
   tile
   topk
   trace
   transpose
   true_divide
   tril
   triu
   trunc
   unique
   unravel_index
   var
   vdot
   vstack
   where
   zeros
   zeros_like

Array
~~~~~

.. autosummary::
   :toctree: generated/

   Array
   Array.all
   Array.any
   Array.argmax
   Array.argmin
   Array.argtopk
   Array.astype
   Array.blocks
   Array.choose
   Array.chunks
   Array.chunksize
   Array.clip
   Array.compute
   Array.compute_chunk_sizes
   Array.conj
   Array.copy
   Array.cumprod
   Array.cumsum
   Array.dask
   Array.dot
   Array.dtype
   Array.flatten
   Array.imag
   Array.itemsize
   Array.map_blocks
   Array.map_overlap
   Array.max
   Array.mean
   Array.min
   Array.moment
   Array.name
   Array.nbytes
   Array.ndim
   Array.nonzero
   Array.npartitions
   Array.numblocks
   Array.partitions
   Array.persist
   Array.prod
   Array.ravel
   Array.real
   Array.rechunk
   Array.repeat
   Array.reshape
   Array.round
   Array.shape
   Array.size
   Array.squeeze
   Array.std
   Array.store
   Array.sum
   Array.swapaxes
   Array.to_dask_dataframe
   Array.to_delayed
   Array.to_hdf5
   Array.to_svg
   Array.to_tiledb
   Array.to_zarr
   Array.topk
   Array.trace
   Array.transpose
   Array.var
   Array.view
   Array.vindex
   Array.visualize


Fast Fourier Transforms
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   fft.fft_wrap
   fft.fft
   fft.fft2
   fft.fftn
   fft.ifft
   fft.ifft2
   fft.ifftn
   fft.rfft
   fft.rfft2
   fft.rfftn
   fft.irfft
   fft.irfft2
   fft.irfftn
   fft.hfft
   fft.ihfft
   fft.fftfreq
   fft.rfftfreq
   fft.fftshift
   fft.ifftshift

Linear Algebra
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   linalg.cholesky
   linalg.inv
   linalg.lstsq
   linalg.lu
   linalg.norm
   linalg.qr
   linalg.solve
   linalg.solve_triangular
   linalg.svd
   linalg.svd_compressed
   linalg.sfqr
   linalg.tsqr

Masked Arrays
~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   ma.average
   ma.filled
   ma.fix_invalid
   ma.getdata
   ma.getmaskarray
   ma.masked_array
   ma.masked_equal
   ma.masked_greater
   ma.masked_greater_equal
   ma.masked_inside
   ma.masked_invalid
   ma.masked_less
   ma.masked_less_equal
   ma.masked_not_equal
   ma.masked_outside
   ma.masked_values
   ma.masked_where
   ma.set_fill_value

Random
~~~~~~

.. autosummary::
   :toctree: generated/

   random.beta
   random.binomial
   random.chisquare
   random.choice
   random.exponential
   random.f
   random.gamma
   random.geometric
   random.gumbel
   random.hypergeometric
   random.laplace
   random.logistic
   random.lognormal
   random.logseries
   random.negative_binomial
   random.noncentral_chisquare
   random.noncentral_f
   random.normal
   random.pareto
   random.permutation
   random.poisson
   random.power
   random.randint
   random.random
   random.random_sample
   random.rayleigh
   random.standard_cauchy
   random.standard_exponential
   random.standard_gamma
   random.standard_normal
   random.standard_t
   random.triangular
   random.uniform
   random.vonmises
   random.wald
   random.weibull
   random.zipf

Stats
~~~~~

.. autosummary::
   :toctree: generated/

   stats.ttest_ind
   stats.ttest_1samp
   stats.ttest_rel
   stats.chisquare
   stats.power_divergence
   stats.skew
   stats.skewtest
   stats.kurtosis
   stats.kurtosistest
   stats.normaltest
   stats.f_oneway
   stats.moment

Image Support
~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   image.imread

Slightly Overlapping Computations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   overlap.overlap
   overlap.map_overlap
   lib.stride_tricks.sliding_window_view
   overlap.trim_internal
   overlap.trim_overlap


Create and Store Arrays
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   from_array
   from_delayed
   from_npy_stack
   from_zarr
   from_tiledb
   store
   to_hdf5
   to_zarr
   to_npy_stack
   to_tiledb

Generalized Ufuncs
~~~~~~~~~~~~~~~~~~

.. currentmodule:: dask.array.gufunc

.. autosummary::
   :toctree: generated/

   apply_gufunc
   as_gufunc
   gufunc


Internal functions
~~~~~~~~~~~~~~~~~~

.. currentmodule:: dask.array.core

.. autosummary::
   :toctree: generated/

   blockwise
   normalize_chunks
   unify_chunks
Python API (advanced)
=====================

.. currentmodule:: distributed

In some rare cases, experts may want to create ``Scheduler``, ``Worker``, and
``Nanny``  objects explicitly in Python.  This is often necessary when making
tools to automatically deploy Dask in custom settings.

It is more common to create a :doc:`Local cluster with Client() on a single
machine <deploying-python>` or use the :doc:`Command Line Interface (CLI) <deploying-cli>`.
New readers are recommended to start there.

If you do want to start Scheduler and Worker objects yourself you should be a
little familiar with ``async``/``await`` style Python syntax.  These objects
are awaitable and are commonly used within ``async with`` context managers.
Here are a few examples to show a few ways to start and finish things.

Full Example
------------

.. autosummary::
   Scheduler
   Worker
   Client

We first start with a comprehensive example of setting up a Scheduler, two Workers,
and one Client in the same event loop, running a simple computation, and then
cleaning everything up.

.. code-block:: python

   import asyncio
   from dask.distributed import Scheduler, Worker, Client

   async def f():
       async with Scheduler() as s:
           async with Worker(s.address) as w1, Worker(s.address) as w2:
               async with Client(s.address, asynchronous=True) as client:
                   future = client.submit(lambda x: x + 1, 10)
                   result = await future
                   print(result)

   asyncio.get_event_loop().run_until_complete(f())

Now we look at simpler examples that build up to this case.

Scheduler
---------

.. autosummary::
   Scheduler

We create scheduler by creating a ``Scheduler()`` object, and then ``await``
that object to wait for it to start up.  We can then wait on the ``.finished``
method to wait until it closes.  In the meantime the scheduler will be active
managing the cluster..

.. code-block:: python

   import asyncio
   from dask.distributed import Scheduler, Worker

   async def f():
       s = Scheduler()        # scheduler created, but not yet running
       s = await s            # the scheduler is running
       await s.finished()     # wait until the scheduler closes

   asyncio.get_event_loop().run_until_complete(f())

This program will run forever, or until some external process connects to the
scheduler and tells it to stop.  If you want to close things yourself you can
close any ``Scheduler``, ``Worker``, ``Nanny``, or ``Client`` class by awaiting
the ``.close`` method:

.. code-block:: python

   await s.close()


Worker
------

.. autosummary::
   Worker

The worker follows the same API.
The only difference is that the worker needs to know the address of the
scheduler.

.. code-block:: python

   import asyncio
   from dask.distributed import Scheduler, Worker

   async def f(scheduler_address):
       w = await Worker(scheduler_address)
       await w.finished()

   asyncio.get_event_loop().run_until_complete(f("tcp://127.0.0.1:8786"))


Start many in one event loop
----------------------------

.. autosummary::
   Scheduler
   Worker

We can run as many of these objects as we like in the same event loop.

.. code-block:: python

   import asyncio
   from dask.distributed import Scheduler, Worker

   async def f():
       s = await Scheduler()
       w = await Worker(s.address)
       await w.finished()
       await s.finished()

   asyncio.get_event_loop().run_until_complete(f())


Use Context Managers
--------------------

We can also use ``async with`` context managers to make sure that we clean up
properly.  Here is the same example as from above:

.. code-block:: python

   import asyncio
   from dask.distributed import Scheduler, Worker

   async def f():
       async with Scheduler() as s:
           async with Worker(s.address) as w:
               await w.finished()
               await s.finished()

   asyncio.get_event_loop().run_until_complete(f())

Alternatively, in the example below we also include a ``Client``, run a small
computation, and then allow things to clean up after that computation..

.. code-block:: python

   import asyncio
   from dask.distributed import Scheduler, Worker, Client

   async def f():
       async with Scheduler() as s:
           async with Worker(s.address) as w1, Worker(s.address) as w2:
               async with Client(s.address, asynchronous=True) as client:
                   future = client.submit(lambda x: x + 1, 10)
                   result = await future
                   print(result)

   asyncio.get_event_loop().run_until_complete(f())

This is equivalent to creating and ``awaiting`` each server, and then calling
``.close`` on each as we leave the context.
In this example we don't wait on ``s.finished()``, so this will terminate
relatively quickly.  You could have called ``await s.finished()`` though if you
wanted this to run forever.

Nanny
-----

.. autosummary::
   Nanny

Alternatively, we can replace ``Worker`` with ``Nanny`` if we want your workers
to be managed in a separate process.  The ``Nanny`` constructor follows the
same API. This allows workers to restart themselves in case of failure. Also,
it provides some additional monitoring, and is useful when coordinating many
workers that should live in different processes in order to avoid the GIL_.

.. code-block:: python

   # w = await Worker(s.address)
   w = await Nanny(s.address)

.. _GIL: https://docs.python.org/3/glossary.html#term-gil


API
---

These classes have a variety of keyword arguments that you can use to control
their behavior.  See the API documentation below for more information.

Scheduler
~~~~~~~~~
.. autoclass:: Scheduler
   :members:

Worker
~~~~~~

.. autoclass:: Worker
   :members:

Nanny
~~~~~

.. autoclass:: Nanny
   :members:
.. _dataframe.design:

Internal Design
===============

Dask DataFrames coordinate many Pandas DataFrames/Series arranged along an
index.  We define a Dask DataFrame object with the following components:

- A Dask graph with a special set of keys designating partitions, such as
  ``('x', 0), ('x', 1), ...``
- A name to identify which keys in the Dask graph refer to this DataFrame, such
  as ``'x'``
- An empty Pandas object containing appropriate metadata (e.g.  column names,
  dtypes, etc.)
- A sequence of partition boundaries along the index called ``divisions``

Metadata
--------

Many DataFrame operations rely on knowing the name and dtype of columns.  To
keep track of this information, all Dask DataFrame objects have a ``_meta``
attribute which contains an empty Pandas object with the same dtypes and names.
For example:

.. code-block:: python

   >>> df = pd.DataFrame({'a': [1, 2, 3], 'b': ['x', 'y', 'z']})
   >>> ddf = dd.from_pandas(df, npartitions=2)
   >>> ddf._meta
   Empty DataFrame
   Columns: [a, b]
   Index: []
   >>> ddf._meta.dtypes
   a     int64
   b    object
   dtype: object

Internally, Dask DataFrame does its best to propagate this information
through all operations, so most of the time a user shouldn't have to worry
about this.  Usually this is done by evaluating the operation on a small sample
of fake data, which can be found on the ``_meta_nonempty`` attribute:

.. code-block:: python

   >>> ddf._meta_nonempty
      a    b
   0  1  foo
   1  1  foo

Sometimes this operation may fail in user defined functions (e.g. when using
``DataFrame.apply``), or may be prohibitively expensive.  For these cases, many
functions support an optional ``meta`` keyword, which allows specifying the
metadata directly, avoiding the inference step.  For convenience, this supports
several options:

1. A Pandas object with appropriate dtypes and names.  If not empty, an empty
   slice will be taken:

.. code-block:: python

  >>> ddf.map_partitions(foo, meta=pd.DataFrame({'a': [1], 'b': [2]}))

2. A description of the appropriate names and dtypes.  This can take several forms:

    * A ``dict`` of ``{name: dtype}`` or an iterable of ``(name, dtype)``
      specifies a DataFrame. Note that order is important: the order of the
      names in ``meta`` should match the order of the columns
    * A tuple of ``(name, dtype)`` specifies a series
    * A dtype object or string (e.g. ``'f8'``) specifies a scalar

This keyword is available on all functions/methods that take user provided
callables (e.g. ``DataFrame.map_partitions``, ``DataFrame.apply``, etc...), as
well as many creation functions (e.g. ``dd.from_delayed``).


.. _dataframe-design-partitions:

Partitions
----------

Internally, a Dask DataFrame is split into many partitions, where each partition
is one Pandas DataFrame.  These DataFrames are split vertically along the
index.  When our index is sorted and we know the values of the divisions of our
partitions, then we can be clever and efficient with expensive algorithms (e.g.
groupby's, joins, etc...).

For example, if we have a time-series index, then our partitions might be
divided by month: all of January will live in one partition while all of
February will live in the next.  In these cases, operations like ``loc``,
``groupby``, and ``join/merge`` along the index can be *much* more efficient
than would otherwise be possible in parallel.  You can view the number of
partitions and divisions of your DataFrame with the following fields:

.. code-block:: python

   >>> df.npartitions
   4
   >>> df.divisions
   ['2015-01-01', '2015-02-01', '2015-03-01', '2015-04-01', '2015-04-31']

Divisions includes the minimum value of every partition's index and the maximum
value of the last partition's index.  In the example above, if the user searches
for a specific datetime range, then we know which partitions we need to inspect
and which we can drop:

.. code-block:: python

   >>> df.loc['2015-01-20': '2015-02-10']  # Must inspect first two partitions

Often we do not have such information about our partitions.  When reading CSV
files, for example, we do not know, without extra user input, how the data is
divided.  In this case ``.divisions`` will be all ``None``:

.. code-block:: python

   >>> df.divisions
   [None, None, None, None, None]

In these cases, any operation that requires a cleanly partitioned DataFrame with
known divisions will have to perform a sort.  This can generally achieved by
calling ``df.set_index(...)``.


.. _dataframe-design-groupby:

Groupby
-------

By default, groupby methods return an object with only 1 partition. This is to
optimize performance, and assumes the groupby reduction returns an object that
is small enough to fit into memory. If your returned object is larger than this,
you can increase the number of output partitions using the `split_out` argument.

.. code-block:: python

   result = df.groupby('id').value.mean()
   result.npartitions  # returns 1

   result = df.groupby('id').value.mean(split_out=8)
   result.npartitions  # returns 8
Stack, Concatenate, and Block
=============================

Often we have many arrays stored on disk that we want to stack together and
think of as one large array.  This is common with geospatial data in which we
might have many HDF5/NetCDF files on disk, one for every day, but we want to do
operations that span multiple days.

To solve this problem, we use the functions ``da.stack``, ``da.concatenate``,
and ``da.block``.

Stack
-----

We stack many existing Dask arrays into a new array, creating a new dimension
as we go.

.. code-block:: python

   >>> import dask.array as da

   >>> arr0 = da.from_array(np.zeros((3, 4)), chunks=(1, 2))
   >>> arr1 = da.from_array(np.ones((3, 4)), chunks=(1, 2))

   >>> data = [arr0, arr1]

   >>> x = da.stack(data, axis=0)
   >>> x.shape
   (2, 3, 4)

   >>> da.stack(data, axis=1).shape
   (3, 2, 4)

   >>> da.stack(data, axis=-1).shape
   (3, 4, 2)

This creates a new dimension with length equal to the number of slices

Concatenate
-----------

We concatenate existing arrays into a new array, extending them along an
existing dimension

.. code-block:: python

   >>> import dask.array as da
   >>> import numpy as np

   >>> arr0 = da.from_array(np.zeros((3, 4)), chunks=(1, 2))
   >>> arr1 = da.from_array(np.ones((3, 4)), chunks=(1, 2))

   >>> data = [arr0, arr1]

   >>> x = da.concatenate(data, axis=0)
   >>> x.shape
   (6, 4)

   >>> da.concatenate(data, axis=1).shape
   (3, 8)

Block
-----

We can handle a larger variety of cases with ``da.block`` as it allows
concatenation to be applied over multiple dimensions at once.  This is useful if
your chunks tile a space, for example if small squares tile a larger 2-D plane.

.. code-block:: python

   >>> import dask.array as da
   >>> import numpy as np

   >>> arr0 = da.from_array(np.zeros((3, 4)), chunks=(1, 2))
   >>> arr1 = da.from_array(np.ones((3, 4)), chunks=(1, 2))

   >>> data = [
   ...     [arr0, arr1],
   ...     [arr1, arr0]
   ... ]

   >>> x = da.block(data)
   >>> x.shape
   (6, 8)
Best Practices
==============

It is easy to get started with Dask delayed, but using it *well* does require
some experience.  This page contains suggestions for best practices, and
includes solutions to common problems.


Call delayed on the function, not the result
--------------------------------------------

Dask delayed operates on functions like ``dask.delayed(f)(x, y)``, not on their results like ``dask.delayed(f(x, y))``.  When you do the latter, Python first calculates ``f(x, y)`` before Dask has a chance to step in.

+------------------------------------------------+--------------------------------------------------------------+
| **Don't**                                      | **Do**                                                       |
+------------------------------------------------+--------------------------------------------------------------+
| .. code-block:: python                         | .. code-block:: python                                       |
|                                                |                                                              |
|    # This executes immediately                 |    # This makes a delayed function, acting lazily            |
|                                                |                                                              |
|    dask.delayed(f(x, y))                       |    dask.delayed(f)(x, y)                                     |
|                                                |                                                              |
+------------------------------------------------+--------------------------------------------------------------+


Compute on lots of computation at once
--------------------------------------

To improve parallelism, you want to include lots of computation in each compute call.
Ideally, you want to make many ``dask.delayed`` calls to define your computation and
then call ``dask.compute`` only at the end.  It is ok to call ``dask.compute``
in the middle of your computation as well, but everything will stop there as
Dask computes those results before moving forward with your code.

+--------------------------------------------------------+-------------------------------------------+
| **Don't**                                              | **Do**                                    |
+--------------------------------------------------------+-------------------------------------------+
| .. code-block:: python                                 | .. code-block:: python                    |
|                                                        |                                           |
|    # Avoid calling compute repeatedly                  |    # Collect many calls for one compute   |
|                                                        |                                           |
|    results = []                                        |    results = []                           |
|    for x in L:                                         |    for x in L:                            |
|        y = dask.delayed(f)(x)                          |        y = dask.delayed(f)(x)             |
|        results.append(y.compute())                     |        results.append(y)                  |
|                                                        |                                           |
|    results                                             |    results = dask.compute(*results)       |
+--------------------------------------------------------+-------------------------------------------+

Calling `y.compute()` within the loop would await the result of the computation every time, and
so inhibit parallelism.

Don't mutate inputs
-------------------

Your functions should not change the inputs directly.

+-----------------------------------------+--------------------------------------+
| **Don't**                               | **Do**                               |
+-----------------------------------------+--------------------------------------+
| .. code-block:: python                  | .. code-block:: python               |
|                                         |                                      |
|    # Mutate inputs in functions         |    # Return new values or copies     |
|                                         |                                      |
|    @dask.delayed                        |    @dask.delayed                     |
|    def f(x):                            |    def f(x):                         |
|        x += 1                           |        x = x + 1                     |
|        return x                         |        return x                      |
+-----------------------------------------+--------------------------------------+

If you need to use a mutable operation, then make a copy within your function first:

.. code-block:: python

   @dask.delayed
   def f(x):
       x = copy(x)
       x += 1
       return x


Avoid global state
------------------

Ideally, your operations shouldn't rely on global state.  Using global state
*might* work if you only use threads, but when you move to multiprocessing or
distributed computing then you will likely encounter confusing errors.

+-------------------------------------------+
| **Don't**                                 |
+-------------------------------------------+
| .. code-block:: python                    |
|                                           |
|   L = []                                  |
|                                           |
|   # This references global variable L     |
|                                           |
|   @dask.delayed                           |
|   def f(x):                               |
|       L.append(x)                         |
|                                           |
+-------------------------------------------+



Don't rely on side effects
--------------------------

Delayed functions only do something if they are computed.  You will always need
to pass the output to something that eventually calls compute.

+--------------------------------+-----------------------------------------+
| **Don't**                      | **Do**                                  |
+--------------------------------+-----------------------------------------+
| .. code-block:: python         | .. code-block:: python                  |
|                                |                                         |
|    # Forget to call compute    |    # Ensure delayed tasks are computed  |
|                                |                                         |
|    dask.delayed(f)(1, 2, 3)    |    x = dask.delayed(f)(1, 2, 3)         |
|                                |    ...                                  |
|    ...                         |    dask.compute(x, ...)                 |
+--------------------------------+-----------------------------------------+

In the first case here, nothing happens, because ``compute()`` is never called.

Break up computations into many pieces
--------------------------------------

Every ``dask.delayed`` function call is a single operation from Dask's perspective.
You achieve parallelism by having many delayed calls, not by using only a
single one: Dask will not look inside a function decorated with ``@dask.delayed``
and parallelize that code internally.  To accomplish that, it needs your help to
find good places to break up a computation.

+------------------------------------+--------------------------------------+
| **Don't**                          | **Do**                               |
+------------------------------------+--------------------------------------+
| .. code-block:: python             | .. code-block:: python               |
|                                    |                                      |
|    # One giant task                |    # Break up into many tasks        |
|                                    |                                      |
|                                    |    @dask.delayed                     |
|    def load(filename):             |    def load(filename):               |
|        ...                         |        ...                           |
|                                    |                                      |
|                                    |    @dask.delayed                     |
|    def process(data):              |    def process(data):                |
|        ...                         |        ...                           |
|                                    |                                      |
|                                    |    @dask.delayed                     |
|    def save(data):                 |    def save(data):                   |
|        ...                         |        ...                           |
|                                    |                                      |
|    @dask.delayed                   |                                      |
|    def f(filenames):               |    def f(filenames):                 |
|        results = []                |        results = []                  |
|        for filename in filenames:  |        for filename in filenames:    |
|            data = load(filename)   |            data = load(filename)     |
|            data = process(data)    |            data = process(data)      |
|            result = save(data)     |            result = save(data)       |
|                                    |                                      |
|        return results              |        return results                |
|                                    |                                      |
|    dask.compute(f(filenames))      |    dask.compute(f(filenames))        |
+------------------------------------+--------------------------------------+

The first version only has one delayed task, and so cannot parallelize.

Avoid too many tasks
--------------------

Every delayed task has an overhead of a few hundred microseconds.  Usually this
is ok, but it can become a problem if you apply ``dask.delayed`` too finely.  In
this case, it's often best to break up your many tasks into batches or use one
of the Dask collections to help you.

+------------------------------------+-------------------------------------------------------------+
| **Don't**                          | **Do**                                                      |
+------------------------------------+-------------------------------------------------------------+
| .. code-block:: python             | .. code-block:: python                                      |
|                                    |                                                             |
|    # Too many tasks                |    # Use collections                                        |
|                                    |                                                             |
|    results = []                    |    import dask.bag as db                                    |
|    for x in range(10000000):       |    b = db.from_sequence(range(10000000), npartitions=1000)  |
|        y = dask.delayed(f)(x)      |    b = b.map(f)                                             |
|        results.append(y)           |    ...                                                      |
|                                    |                                                             |
+------------------------------------+-------------------------------------------------------------+

Here we use ``dask.bag`` to automatically batch applying our function. We could also have constructed
our own batching as follows

.. code-block:: python

   def batch(seq):
       sub_results = []
       for x in seq:
           sub_results.append(f(x))
       return sub_results

    batches = []
    for i in range(0, 10000000, 10000):
        result_batch = dask.delayed(batch)(range(i, i + 10000))
        batches.append(result_batch)


Here we construct batches where each delayed function call computes for many data points from
the original input.

Avoid calling delayed within delayed functions
----------------------------------------------

Often, if you are new to using Dask delayed, you place ``dask.delayed`` calls
everywhere and hope for the best.  While this may actually work, it's usually
slow and results in hard-to-understand solutions.

Usually you never call ``dask.delayed`` within ``dask.delayed`` functions.

+----------------------------------------+--------------------------------------+
| **Don't**                              | **Do**                               |
+----------------------------------------+--------------------------------------+
| .. code-block:: python                 | .. code-block:: python               |
|                                        |                                      |
|    # Delayed function calls delayed    |    # Normal function calls delayed   |
|                                        |                                      |
|    @dask.delayed                       |                                      |
|    def process_all(L):                 |    def process_all(L):               |
|        result = []                     |        result = []                   |
|        for x in L:                     |        for x in L:                   |
|            y = dask.delayed(f)(x)      |            y = dask.delayed(f)(x)    |
|            result.append(y)            |            result.append(y)          |
|        return result                   |        return result                 |
+----------------------------------------+--------------------------------------+

Because the normal function only does delayed work it is very fast and so
there is no reason to delay it.

Don't call dask.delayed on other Dask collections
-------------------------------------------------

When you place a Dask array or Dask DataFrame into a delayed call, that function
will receive the NumPy or Pandas equivalent.  Beware that if your array is
large, then this might crash your workers.

Instead, it's more common to use methods like ``da.map_blocks``

+--------------------------------------------------+---------------------------------------------+
| **Don't**                                        | **Do**                                      |
+--------------------------------------------------+---------------------------------------------+
| .. code-block:: python                           | .. code-block:: python                      |
|                                                  |                                             |
|    # Call delayed functions on Dask collections  |    # Use mapping methods if applicable      |
|                                                  |                                             |
|    import dask.dataframe as dd                   |    import dask.dataframe as dd              |
|    df = dd.read_csv('/path/to/*.csv')            |    df = dd.read_csv('/path/to/*.csv')       |
|                                                  |                                             |
|    dask.delayed(train)(df)                       |    df.map_partitions(train)                 |
+--------------------------------------------------+---------------------------------------------+

Alternatively, if the procedure doesn't fit into a mapping, you can always
turn your arrays or dataframes into *many* delayed
objects, for example

.. code-block:: python

    partitions = df.to_delayed()
    delayed_values = [dask.delayed(train)(part)
                      for part in partitions]

However, if you don't mind turning your Dask array/DataFrame into a single
chunk, then this is ok.

.. code-block:: python

   dask.delayed(train)(..., y=df.sum())


Avoid repeatedly putting large inputs into delayed calls
--------------------------------------------------------

Every time you pass a concrete result (anything that isn't delayed) Dask will
hash it by default to give it a name.  This is fairly fast (around 500 MB/s)
but can be slow if you do it over and over again.  Instead, it is better to
delay your data as well.

This is especially important when using a distributed cluster to avoid sending
your data separately for each function call.

+------------------------------------------+---------------------------------------------------------+
| **Don't**                                | **Do**                                                  |
+------------------------------------------+---------------------------------------------------------+
| .. code-block:: python                   | .. code-block:: python                                  |
|                                          |                                                         |
|    x = np.array(...)  # some large array |    x = np.array(...)    # some large array              |
|                                          |    x = dask.delayed(x)  # delay the data once           |
|    results = [dask.delayed(train)(x, i)  |    results = [dask.delayed(train)(x, i)                 |
|               for i in range(1000)]      |               for i in range(1000)]                     |
+------------------------------------------+---------------------------------------------------------+


Every call to ``dask.delayed(train)(x, ...)`` has to hash the NumPy array ``x``, which slows things down.
.. _scheduling-policy:

Scheduling in Depth
===================

*Note: this technical document is not optimized for user readability.*

The default shared memory scheduler used by most dask collections lives in
``dask/local.py``. This scheduler dynamically schedules tasks to new
workers as they become available.  It operates in a shared memory environment
without consideration to data locality, all workers have access to all data
equally.

We find that our workloads are best served by trying to minimize the memory
footprint.  This document talks about our policies to accomplish this in our
scheduling budget of one millisecond per task, irrespective of the number of
tasks.

Generally we are faced with the following situation:  A worker arrives with a
newly completed task.  We update our data structures of execution state and
have to provide a new task for that worker.  In general there are very many
available tasks, which should we give to the worker?

*Q: Which of our available tasks should we give to the newly ready worker?*

This question is simple and local and yet strongly impacts the performance of
our algorithm.  We want to choose a task that lets us free memory now and in
the future.  We need a clever and cheap way to break a tie between the set of
available tasks.

At this stage we choose the policy of "last in, first out."  That is we choose
the task that was most recently made available, quite possibly by the worker
that just returned to us.  This encourages the general theme of finishing
things before starting new things.

We implement this with a stack.  When a worker arrives with its finished task
we figure out what new tasks we can now compute with the new data and put those
on top of the stack if any exist.  We pop an item off of the top of the stack
and deliver that to the waiting worker.

And yet if the newly completed task makes ready multiple newly ready tasks in
which order should we place them on the stack?  This is yet another opportunity
for a tie breaker.  This is particularly important at *the beginning* of
execution where we typically add a large number of leaf tasks onto the stack.
Our choice in this tie breaker also strongly affects performance in many cases.

We want to encourage depth first behavior where, if our computation is composed
of something like many trees we want to fully explore one subtree before moving
on to the next.  This encourages our workers to complete blocks/subtrees of our
graph before moving on to new blocks/subtrees.

And so to encourage this "depth first behavior" we do a depth first search and
number all nodes according to their number in the depth first search (DFS)
traversal.  We use this number to break ties when adding tasks on to the stack.
Please note that while we spoke of optimizing the many-distinct-subtree case
above this choice is entirely local and applies quite generally beyond this
case.  Anything that behaves even remotely like the many-distinct-subtree case
will benefit accordingly, and this case is quite common in normal workloads.

And yet we have glossed over another tie breaker. Performing the depth
first search, when we arrive at a node with many children we can choose the
order in which to traverse the children.  We resolve this tie breaker by
selecting those children whose result is depended upon by the most nodes.  This
dependence can be either direct for those nodes that take that data as input or
indirect for any ancestor node in the graph.  This emphasizing traversing first
those nodes that are parts of critical paths having long vertical chains that
rest on top of this node's result, and nodes whose data is depended upon by
many nodes in the future.  We choose to dive down into these subtrees first in
our depth first search so that future computations don't get stuck waiting for
them to complete.

And so we have three tie breakers

1.  Q:  Which of these available tasks should I run?

    A:  Last in, first out
2.  Q:  Which of these tasks should I put on the stack first?

    A:  Do a depth first search before the computation, use that ordering.
3.  Q:  When performing the depth first search how should I choose between
    children?

    A:  Choose those children on whom the most data depends

We have found common workflow types that require each of these decisions.  We
have not yet run into a commonly occurring graph type in data analysis that is
not well handled by these heuristics for the purposes of minimizing memory use.
Futures
=======

Dask supports a real-time task framework that extends Python's
`concurrent.futures <https://docs.python.org/3/library/concurrent.futures.html>`_
interface.  This interface is good for arbitrary task scheduling like
:doc:`dask.delayed <delayed>`, but is immediate rather than lazy, which
provides some more flexibility in situations where the computations may evolve
over time.

These features depend on the second generation task scheduler found in
`dask.distributed <https://distributed.dask.org/en/latest>`_ (which,
despite its name, runs very well on a single machine).

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/07EiCpdhtDE"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

.. currentmodule:: distributed

Examples
--------

Visit https://examples.dask.org/futures.html to see and run examples
using futures with Dask.

Start Dask Client
-----------------

You must start a ``Client`` to use the futures interface.  This tracks state
among the various worker processes or threads:

.. code-block:: python

   from dask.distributed import Client

   client = Client()  # start local workers as processes
   # or
   client = Client(processes=False)  # start local workers as threads

If you have `Bokeh <https://docs.bokeh.org>`_ installed, then this starts up a
diagnostic dashboard at ``http://localhost:8787`` .

Submit Tasks
------------

.. autosummary::
   Client.submit
   Client.map
   Future.result

You can submit individual tasks using the ``submit`` method:

.. code-block:: python

   def inc(x):
       return x + 1

   def add(x, y):
       return x + y

   a = client.submit(inc, 10)  # calls inc(10) in background thread or process
   b = client.submit(inc, 20)  # calls inc(20) in background thread or process

The ``submit`` function returns a ``Future``, which refers to a remote result.  This result may
not yet be completed:

.. code-block:: python

   >>> a
   <Future: status: pending, key: inc-b8aaf26b99466a7a1980efa1ade6701d>

Eventually it will complete.  The result stays in the remote
thread/process/worker until you ask for it back explicitly:

.. code-block:: python

   >>> a
   <Future: status: finished, type: int, key: inc-b8aaf26b99466a7a1980efa1ade6701d>

   >>> a.result()  # blocks until task completes and data arrives
   11

You can pass futures as inputs to submit.  Dask automatically handles dependency
tracking; once all input futures have completed, they will be moved onto a
single worker (if necessary), and then the computation that depends on them
will be started.  You do not need to wait for inputs to finish before
submitting a new task; Dask will handle this automatically:

.. code-block:: python

   c = client.submit(add, a, b)  # calls add on the results of a and b

Similar to Python's ``map``, you can use ``Client.map`` to call the same
function and many inputs:

.. code-block:: python

   futures = client.map(inc, range(1000))

However, note that each task comes with about 1ms of overhead.  If you want to
map a function over a large number of inputs, then you might consider
:doc:`dask.bag <bag>` or :doc:`dask.dataframe <dataframe>` instead.

.. note: See `this page <https://docs.dask.org/en/latest/graphs.html>`_ for
   restrictions on what functions you use with Dask.

Move Data
---------

.. autosummary::
   Future.result
   Client.gather
   Client.scatter

Given any future, you can call the ``.result`` method to gather the result.
This will block until the future is done computing and then transfer the result
back to your local process if necessary:

.. code-block:: python

   >>> c.result()
   32

You can gather many results concurrently using the ``Client.gather`` method.
This can be more efficient than calling ``.result()`` on each future
sequentially:

.. code-block:: python

   >>> # results = [future.result() for future in futures]
   >>> results = client.gather(futures)  # this can be faster

If you have important local data that you want to include in your computation,
you can either include it as a normal input to a submit or map call:

.. code-block:: python

   >>> df = pd.read_csv('training-data.csv')
   >>> future = client.submit(my_function, df)

Or you can ``scatter`` it explicitly.  Scattering moves your data to a worker
and returns a future pointing to that data:

.. code-block:: python

   >>> remote_df = client.scatter(df)
   >>> remote_df
   <Future: status: finished, type: DataFrame, key: bbd0ca93589c56ea14af49cba470006e>

   >>> future = client.submit(my_function, remote_df)

Both of these accomplish the same result, but using scatter can sometimes be
faster.  This is especially true if you use processes or distributed workers
(where data transfer is necessary) and you want to use ``df`` in many
computations.  Scattering the data beforehand avoids excessive data movement.

Calling scatter on a list scatters all elements individually.  Dask will spread
these elements evenly throughout workers in a round-robin fashion:

.. code-block:: python

   >>> client.scatter([1, 2, 3])
   [<Future: status: finished, type: int, key: c0a8a20f903a4915b94db8de3ea63195>,
    <Future: status: finished, type: int, key: 58e78e1b34eb49a68c65b54815d1b158>,
    <Future: status: finished, type: int, key: d3395e15f605bc35ab1bac6341a285e2>]

References, Cancellation, and Exceptions
----------------------------------------

.. autosummary::
   Future.cancel
   Future.exception
   Future.traceback
   Client.cancel

Dask will only compute and hold onto results for which there are active
futures.  In this way, your local variables define what is active in Dask.  When
a future is garbage collected by your local Python session, Dask will feel free
to delete that data or stop ongoing computations that were trying to produce
it:

.. code-block:: python

   >>> del future  # deletes remote data once future is garbage collected

You can also explicitly cancel a task using the ``Future.cancel`` or
``Client.cancel`` methods:

.. code-block:: python

   >>> future.cancel()  # deletes data even if other futures point to it

If a future fails, then Dask will raise the remote exceptions and tracebacks if
you try to get the result:

.. code-block:: python

   def div(x, y):
       return x / y

   >>> a = client.submit(div, 1, 0)  # 1 / 0 raises a ZeroDivisionError
   >>> a
   <Future: status: error, key: div-3601743182196fb56339e584a2bf1039>

   >>> a.result()
         1 def div(x, y):
   ----> 2     return x / y

   ZeroDivisionError: division by zero

All futures that depend on an erred future also err with the same exception:

.. code-block:: python

   >>> b = client.submit(inc, a)
   >>> b
   <Future: status: error, key: inc-15e2e4450a0227fa38ede4d6b1a952db>

You can collect the exception or traceback explicitly with the
``Future.exception`` or ``Future.traceback`` methods.


Waiting on Futures
------------------

.. autosummary::
   as_completed
   wait

You can wait on a future or collection of futures using the ``wait`` function:

.. code-block:: python

   from dask.distributed import wait

   >>> wait(futures)

This blocks until all futures are finished or have erred.

You can also iterate over the futures as they complete using the
``as_completed`` function:

.. code-block:: python

   from dask.distributed import as_completed

   futures = client.map(score, x_values)

   best = -1
   for future in as_completed(futures):
      y = future.result()
      if y > best:
          best = y

For greater efficiency, you can also ask ``as_completed`` to gather the results
in the background:

.. code-block:: python

   for future, result in as_completed(futures, with_results=True):
       # y = future.result()  # don't need this
      ...

Or collect all futures in batches that had arrived since the last iteration:

.. code-block:: python

   for batch in as_completed(futures, with_results=True).batches():
      for future, result in batch:
          ...

Additionally, for iterative algorithms, you can add more futures into the
``as_completed`` iterator *during* iteration:

.. code-block:: python

   seq = as_completed(futures)

   for future in seq:
       y = future.result()
       if condition(y):
           new_future = client.submit(...)
           seq.add(new_future)  # add back into the loop

or use ``seq.update(futures)`` to add multiple futures at once.


Fire and Forget
---------------

.. autosummary::
   fire_and_forget

Sometimes we don't care about gathering the result of a task, and only care
about side effects that it might have like writing a result to a file:

.. code-block:: python

   >>> a = client.submit(load, filename)
   >>> b = client.submit(process, a)
   >>> c = client.submit(write, b, out_filename)

As noted above, Dask will stop work that doesn't have any active futures.  It
thinks that because no one has a pointer to this data that no one cares.  You
can tell Dask to compute a task anyway, even if there are no active futures,
using the ``fire_and_forget`` function:

.. code-block:: python

   from dask.distributed import fire_and_forget

   >>> fire_and_forget(c)

This is particularly useful when a future may go out of scope, for example, as
part of a function:

.. code-block:: python

    def process(filename):
        out_filename = 'out-' + filename
        a = client.submit(load, filename)
        b = client.submit(process, a)
        c = client.submit(write, b, out_filename)
        fire_and_forget(c)
        return  # here we lose the reference to c, but that's now ok

    for filename in filenames:
        process(filename)


Submit Tasks from Tasks
-----------------------

.. autosummary::
   get_client
   rejoin
   secede

*This is an advanced feature and is rarely necessary in the common case.*

Tasks can launch other tasks by getting their own client.  This enables complex
and highly dynamic workloads:

.. code-block:: python

   from dask.distributed import get_client

   def my_function(x):
       ...

       # Get locally created client
       client = get_client()

       # Do normal client operations, asking cluster for computation
       a = client.submit(...)
       b = client.submit(...)
       a, b = client.gather([a, b])

       return a + b

It also allows you to set up long running tasks that watch other resources like
sockets or physical sensors:

.. code-block:: python

   def monitor(device):
      client = get_client()
      while True:
          data = device.read_data()
          future = client.submit(process, data)
          fire_and_forget(future)

   for device in devices:
       fire_and_forget(client.submit(monitor))

However, each running task takes up a single thread, and so if you launch many
tasks that launch other tasks, then it is possible to deadlock the system if you
are not careful.  You can call the ``secede`` function from within a task to
have it remove itself from the dedicated thread pool into an administrative
thread that does not take up a slot within the Dask worker:

.. code-block:: python

   from dask.distributed import get_client, secede

   def monitor(device):
      client = get_client()
      secede()  # remove this task from the thread pool
      while True:
          data = device.read_data()
          future = client.submit(process, data)
          fire_and_forget(future)

If you intend to do more work in the same thread after waiting on client work,
you may want to explicitly block until the thread is able to *rejoin* the
thread pool.  This allows some control over the number of threads that are
created and stops too many threads from being active at once, over-saturating
your hardware:

.. code-block:: python

   def f(n):  # assume that this runs as a task
      client = get_client()

      secede()  # secede while we wait for results to come back
      futures = client.map(func, range(n))
      results = client.gather(futures)

      rejoin()  # block until a slot is open in the thread pool
      result = analyze(results)
      return result


Alternatively, you can just use the normal ``compute`` function *within* a
task.  This will automatically call ``secede`` and ``rejoin`` appropriately:

.. code-block:: python

   def f(name, fn):
       df = dd.read_csv(fn)  # note that this is a dask collection
       result = df[df.name == name].count()

       # This calls secede
       # Then runs the computation on the cluster (including this worker)
       # Then blocks on rejoin, and finally delivers the answer
       result = result.compute()

       return result


Coordination Primitives
-----------------------

.. autosummary::
   Queue
   Variable
   Lock
   Event
   Semaphore
   Pub
   Sub

.. note: These are advanced features and are rarely necessary in the common case.

Sometimes situations arise where tasks, workers, or clients need to coordinate
with each other in ways beyond normal task scheduling with futures.  In these
cases Dask provides additional primitives to help in complex situations.

Dask provides distributed versions of coordination primitives like locks, events,
queues, global variables, and pub-sub systems that, where appropriate, match
their in-memory counterparts.  These can be used to control access to external
resources, track progress of ongoing computations, or share data in
side-channels between many workers, clients, and tasks sensibly.

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/Q-Y3BR1u7c0"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

These features are rarely necessary for common use of Dask.  We recommend that
beginning users stick with using the simpler futures found above (like
``Client.submit`` and ``Client.gather``) rather than embracing needlessly
complex techniques.


Queues
~~~~~~

.. autosummary::
   Queue

Dask queues follow the API for the standard Python Queue, but now move futures
or small messages between clients.  Queues serialize sensibly and reconnect
themselves on remote clients if necessary:

.. code-block:: python

   from dask.distributed import Queue

   def load_and_submit(filename):
       data = load(filename)
       client = get_client()
       future = client.submit(process, data)
       queue.put(future)

   client = Client()

   queue = Queue()

   for filename in filenames:
       future = client.submit(load_and_submit, filename)
       fire_and_forget(future)

   while True:
       future = queue.get()
       print(future.result())


Queues can also send small pieces of information, anything that is msgpack
encodable (ints, strings, bools, lists, dicts, etc.).  This can be useful to
send back small scores or administrative messages:

.. code-block:: python

   def func(x):
       try:
          ...
       except Exception as e:
           error_queue.put(str(e))

   error_queue = Queue()

Queues are mediated by the central scheduler, and so they are not ideal for
sending large amounts of data (everything you send will be routed through a
central point).  They are well suited to move around small bits of metadata, or
futures.  These futures may point to much larger pieces of data safely:

.. code-block:: python

   >>> x = ... # my large numpy array

   # Don't do this!
   >>> q.put(x)

   # Do this instead
   >>> future = client.scatter(x)
   >>> q.put(future)

   # Or use futures for metadata
   >>> q.put({'status': 'OK', 'stage=': 1234})

If you're looking to move large amounts of data between workers, then you might
also want to consider the Pub/Sub system described a few sections below.

Global Variables
~~~~~~~~~~~~~~~~

.. autosummary::
   Variable

Variables are like Queues in that they communicate futures and small data
between clients.  However, variables hold only a single value.  You can get or
set that value at any time:

.. code-block:: python

   >>> var = Variable('stopping-criterion')
   >>> var.set(False)

   >>> var.get()
   False

This is often used to signal stopping criteria or current parameters
between clients.

If you want to share large pieces of information, then scatter the data first:

.. code-block:: python

   >>> parameters = np.array(...)
   >>> future = client.scatter(parameters)
   >>> var.set(future)


Locks
~~~~~

.. autosummary::
   Lock

You can also hold onto cluster-wide locks using the ``Lock`` object.
Dask Locks have the same API as normal ``threading.Lock`` objects, except that
they work across the cluster:

.. code-block:: python

       from dask.distributed import Lock
       lock = Lock()

       with lock:
           # access protected resource

You can manage several locks at the same time.  Lock can either be given a
consistent name or you can pass the lock object around itself.

Using a consistent name is convenient when you want to lock some known named resource:

.. code-block:: python

   from dask.distributed import Lock

   def load(fn):
       with Lock('the-production-database'):
           # read data from filename using some sensitive source
           return ...

   futures = client.map(load, filenames)

Passing around a lock works as well and is easier when you want to create short-term
locks for a particular situation:

.. code-block:: python

   from dask.distributed import Lock
   lock = Lock()

   def load(fn, lock=None):
       with lock:
           # read data from filename using some sensitive source
           return ...

   futures = client.map(load, filenames, lock=lock)

This can be useful if you want to control concurrent access to some external
resource like a database or un-thread-safe library.


Events
~~~~~~

.. autosummary::
   Event

Dask Events mimic ``asyncio.Event`` objects, but on a cluster scope.
They hold a single flag which can be set or cleared.
Clients can wait until the event flag is set.
Different from a ``Lock``, every client can set or clear the flag and there
is no "ownership" of an event.

You can use events to e.g. synchronize multiple clients:

.. code-block:: python

   # One one client
   from dask.distributed import Event

   event = Event("my-event-1")
   event.wait()

The call to wait will block until the event is set, e.g. in another client

.. code-block:: python

   # In another client
   from dask.distributed import Event

   event = Event("my-event-1")

   # do some work

   event.set()

Events can be set, cleared and waited on multiple times.
Every waiter referencing the same event name will be notified on event set
(and not only the first one as in the case of a lock):

.. code-block:: python

   from dask.distributed import Event

   def wait_for_event(x):
      event = Event("my-event")

      event.wait()
      # at this point, all function calls
      # are in sync once the event is set

   futures = client.map(wait_for_event, range(10))

   Event("my-event").set()
   client.gather(futures)


Semaphore
~~~~~~~~~

.. autosummary::
   Semaphore

Similar to the single-valued ``Lock`` it is also possible to use a cluster-wide
semaphore to coordinate and limit access to a sensitive resource like a
database.

.. code-block:: python

   from dask.distributed import Semaphore

   sem = Semaphore(max_leases=2, name="database")

   def access_limited(val, sem):
      with sem:
         # Interact with the DB
         return

   futures = client.map(access_limited, range(10), sem=sem)
   client.gather(futures)
   sem.close()


Publish-Subscribe
~~~~~~~~~~~~~~~~~

.. autosummary::
   Pub
   Sub

Dask implements the `Publish Subscribe pattern <https://en.wikipedia.org/wiki/Publish%E2%80%93subscribe_pattern>`_,
providing an additional channel of communication between ongoing tasks.

.. autoclass:: Pub
   :members:

Actors
------

Actors allow workers to manage rapidly changing state without coordinating with
the central scheduler.  This has the advantage of reducing latency
(worker-to-worker roundtrip latency is around 1ms), reducing pressure on the
centralized scheduler (workers can coordinate actors entirely among each other),
and also enabling workflows that require stateful or in-place memory
manipulation.

However, these benefits come at a cost.  The scheduler is unaware of actors and
so they don't benefit from diagnostics, load balancing, or resilience.  Once an
actor is running on a worker it is forever tied to that worker.  If that worker
becomes overburdened or dies, then there is no opportunity to recover the
workload.

*Because Actors avoid the central scheduler they can be high-performing, but not resilient.*

Example: Counter
~~~~~~~~~~~~~~~~

An actor is a class containing both state and methods that is submitted to a
worker:

.. code-block:: python

   class Counter:
       n = 0

       def __init__(self):
           self.n = 0

       def increment(self):
           self.n += 1
           return self.n

   from dask.distributed import Client
   client = Client()

   future = client.submit(Counter, actor=True)
   counter = future.result()

   >>> counter
   <Actor: Counter, key=Counter-afa1cdfb6b4761e616fa2cfab42398c8>

Method calls on this object produce ``ActorFutures``, which are similar to
normal Futures, but interact only with the worker holding the Actor:

.. code-block:: python

   >>> future = counter.increment()
   >>> future
   <ActorFuture>

   >>> future.result()
   1

Attribute access is synchronous and blocking:

.. code-block:: python

   >>> counter.n
   1


Example: Parameter Server
~~~~~~~~~~~~~~~~~~~~~~~~~

This example will perform the following minimization with a parameter server:

.. math::

   \min_{p\in\mathbb{R}^{1000}} \sum_{i=1}^{1000} (p_i - 1)^2

This is a simple minimization that will serve as an illustrative example.

The Dask Actor will serve as the parameter server that will hold the model.
The client will calculate the gradient of the loss function above.

.. code-block:: python

   import numpy as np

   from dask.distributed import Client
   client = Client(processes=False)

   class ParameterServer:
       def __init__(self):
           self.data = dict()

       def put(self, key, value):
           self.data[key] = value

       def get(self, key):
           return self.data[key]

   def train(params, lr=0.1):
       grad = 2 * (params - 1)  # gradient of (params - 1)**2
       new_params = params - lr * grad
       return new_params

   ps_future = client.submit(ParameterServer, actor=True)
   ps = ps_future.result()

   ps.put('parameters', np.random.random(1000))
   for k in range(20):
       params = ps.get('parameters').result()
       new_params = train(params)
       ps.put('parameters', new_params)
       print(new_params.mean())
       # k=0: "0.5988202981316124"
       # k=10: "0.9569236575164062"

This example works, and the loss function is minimized. The (simple) equation
above is minimize, so each :math:`p_i` converges to 1. If desired, this example
could be adapted to machine learning with a more complex function to minimize.

Asynchronous Operation
~~~~~~~~~~~~~~~~~~~~~~

All operations that require talking to the remote worker are awaitable:

.. code-block:: python

   async def f():
       future = client.submit(Counter, actor=True)
       counter = await future  # gather actor object locally

       counter.increment()  # send off a request asynchronously
       await counter.increment()  # or wait until it was received

       n = await counter.n  # attribute access also must be awaited


API
---

**Client**

.. autosummary::
   Client
   Client.cancel
   Client.compute
   Client.gather
   Client.get
   Client.get_dataset
   Client.get_executor
   Client.has_what
   Client.list_datasets
   Client.map
   Client.ncores
   Client.persist
   Client.profile
   Client.publish_dataset
   Client.rebalance
   Client.replicate
   Client.restart
   Client.run
   Client.run_on_scheduler
   Client.scatter
   Client.shutdown
   Client.scheduler_info
   Client.start_ipython_workers
   Client.start_ipython_scheduler
   Client.submit
   Client.unpublish_dataset
   Client.upload_file
   Client.who_has

**Future**

.. autosummary::
   Future
   Future.add_done_callback
   Future.cancel
   Future.cancelled
   Future.done
   Future.exception
   Future.result
   Future.traceback

**Functions**

.. autosummary::
   as_completed
   fire_and_forget
   get_client
   secede
   rejoin
   wait

.. autofunction:: as_completed
.. autofunction:: fire_and_forget
.. autofunction:: get_client
.. autofunction:: secede
.. autofunction:: rejoin
.. autofunction:: wait

.. autoclass:: Client
   :members:

.. autoclass:: Future
   :members:

.. autoclass:: Queue
   :members:

.. autoclass:: Variable
   :members:

.. autoclass:: Lock
   :members:

.. autoclass:: Event
   :members:

.. autoclass:: Pub
   :members:

.. autoclass:: Sub
   :members:
Shuffling for GroupBy and Join
==============================

.. currentmodule:: dask.dataframe

Operations like ``groupby``, ``join``, and ``set_index`` have special
performance considerations that are different from normal Pandas due to the
parallel, larger-than-memory, and distributed nature of Dask DataFrame.

Easy Case
---------

To start off, common groupby operations like
``df.groupby(columns).reduction()`` for known reductions like ``mean, sum, std,
var, count, nunique`` are all quite fast and efficient, even if partitions are
not cleanly divided with known divisions.  This is the common case.

Additionally, if divisions are known, then applying an arbitrary function to
groups is efficient when the grouping columns include the index.

Joins are also quite fast when joining a Dask DataFrame to a Pandas DataFrame
or when joining two Dask DataFrames along their index.  No special
considerations need to be made when operating in these common cases.

So, if you're doing common groupby and join operations, then you can stop reading
this.  Everything will scale nicely.  Fortunately, this is true most of the
time:

.. code-block:: python

   >>> df.groupby(columns).known_reduction()             # Fast and common case
   >>> df.groupby(columns_with_index).apply(user_fn)     # Fast and common case
   >>> dask_df.join(pandas_df, on=column)                # Fast and common case
   >>> lhs.join(rhs)                                     # Fast and common case
   >>> lhs.merge(rhs, on=columns_with_index)             # Fast and common case

Difficult Cases
---------------

In some cases, such as when applying an arbitrary function to groups (when not
grouping on index with known divisions), when joining along non-index columns,
or when explicitly setting an unsorted column to be the index, we may need to
trigger a full dataset shuffle:

.. code-block:: python

   >>> df.groupby(columns_no_index).apply(user_fn)   # Requires shuffle
   >>> lhs.join(rhs, on=columns_no_index)            # Requires shuffle
   >>> df.set_index(column)                          # Requires shuffle

A shuffle is necessary when we need to re-sort our data along a new index.  For
example, if we have banking records that are organized by time and we now want
to organize them by user ID, then we'll need to move a lot of data around.  In
Pandas all of this data fits in memory, so this operation was easy.  Now that we
don't assume that all data fits in memory, we must be a bit more careful.

Re-sorting the data can be avoided by restricting yourself to the easy cases
mentioned above.

Shuffle Methods
---------------

There are currently two strategies to shuffle data depending on whether you are
on a single machine or on a distributed cluster: shuffle on disk and shuffle 
over the network.

Shuffle on Disk
```````````````

When operating on larger-than-memory data on a single machine, we shuffle by
dumping intermediate results to disk.  This is done using the partd_ project
for on-disk shuffles.

.. _partd: https://github.com/dask/partd

Shuffle over the Network
````````````````````````

When operating on a distributed cluster, the Dask workers may not have access to
a shared hard drive.  In this case, we shuffle data by breaking input partitions
into many pieces based on where they will end up and moving these pieces
throughout the network.  This prolific expansion of intermediate partitions
can stress the task scheduler.  To manage for many-partitioned datasets we
sometimes shuffle in stages, causing undue copies but reducing the ``n**2``
effect of shuffling to something closer to ``n log(n)`` with ``log(n)`` copies.

Selecting methods
`````````````````

Dask will use on-disk shuffling by default, but will switch to task-based
distributed shuffling if the default scheduler is set to use a
``dask.distributed.Client``, such as would be the case if the user sets the
Client as default:

.. code-block:: python

    client = Client('scheduler:8786', set_as_default=True)

Alternatively, if you prefer to avoid defaults, you can configure the global
shuffling method by using the ``dask.config.set(shuffle=...)`` command.
This can be done globally:

.. code-block:: python

    dask.config.set(shuffle='tasks')

    df.groupby(...).apply(...)

or as a context manager:

.. code-block:: python

    with dask.config.set(shuffle='tasks'):
        df.groupby(...).apply(...)


In addition, ``set_index`` also accepts a ``shuffle`` keyword argument that
can be used to select either on-disk or task-based shuffling:

.. code-block:: python

    df.set_index(column, shuffle='disk')
    df.set_index(column, shuffle='tasks')


.. _dataframe.groupby.aggregate:

Aggregate
=========

Dask supports Pandas' ``aggregate`` syntax to run multiple reductions on the
same groups.  Common reductions such as ``max``, ``sum``, ``list`` and ``mean`` are 
directly supported:

.. code-block:: python

    >>> df.groupby(columns).aggregate(['sum', 'mean', 'max', 'min', list])

Dask also supports user defined reductions.  To ensure proper performance, the
reduction has to be formulated in terms of three independent steps. The
``chunk`` step is applied to each partition independently and reduces the data
within a partition. The ``aggregate`` combines the within partition results.
The optional ``finalize`` step combines the results returned from the
``aggregate`` step and should return a single final column. For Dask to
recognize the reduction, it has to be passed as an instance of
``dask.dataframe.Aggregation``.

For example, ``sum`` could be implemented as:

.. code-block:: python

    custom_sum = dd.Aggregation('custom_sum', lambda s: s.sum(), lambda s0: s0.sum())
    df.groupby('g').agg(custom_sum)

The name argument should be different from existing reductions to avoid data
corruption.  The arguments to each function are pre-grouped series objects,
similar to ``df.groupby('g')['value']``.

Many reductions can only be implemented with multiple temporaries. To implement
these reductions, the steps should return tuples and expect multiple arguments.
A mean function can be implemented as:

.. code-block:: python

    custom_mean = dd.Aggregation(
        'custom_mean',
        lambda s: (s.count(), s.sum()),
        lambda count, sum: (count.sum(), sum.sum()),
        lambda count, sum: sum / count,
    )
    df.groupby('g').agg(custom_mean)


For example, let's compute the group-wise extent (maximum - minimum)
for a DataFrame.

.. code-block:: python

   >>> df = pd.DataFrame({
   ...   'a': ['a', 'b', 'a', 'a', 'b'],
   ...   'b': [0, 1, 0, 2, 5],
   ... })
   >>> ddf = dd.from_pandas(df, 2)

We define the building blocks to find the maximum and minimum of each chunk, and then
the maximum and minimum over all the chunks. We finalize by taking the difference between
the Series with the maxima and minima

.. code-block:: python

   >>> def chunk(grouped):
   ...     return grouped.max(), grouped.min()

   >>> def agg(chunk_maxes, chunk_mins):
   ...     return chunk_maxes.max(), chunk_mins.min()

   >>> def finalize(maxima, minima):
   ...     return maxima - minima

Finally, we create and use the aggregation

.. code-block:: python

   >>> extent = dd.Aggregation('extent', chunk, agg, finalize=finalize)
   >>> ddf.groupby('a').agg(extent).compute()
      b
   a
   a  2
   b  4

To apply :py:class:`dask.dataframe.groupby.SeriesGroupBy.nunique` to more than one
column you can use:

.. code-block:: python

    >>> df['c'] = [1, 2, 1, 1, 2]
    >>> ddf = dd.from_pandas(df, 2)
    >>> nunique = dd.Aggregation(
    ...     name="nunique",
    ...     chunk=lambda s: s.apply(lambda x: list(set(x))),
    ...     agg=lambda s0: s0.obj.groupby(level=list(range(s0.obj.index.nlevels))).sum(),
    ...     finalize=lambda s1: s1.apply(lambda final: len(set(final))),
    ... )
    >>> ddf.groupby('a').agg({'b':nunique, 'c':nunique})Delayed
=======

.. toctree::
   :maxdepth: 1
   :hidden:

   delayed-collections.rst
   delayed-best-practices.rst

Sometimes problems don't fit into one of the collections like ``dask.array`` or
``dask.dataframe``. In these cases, users can parallelize custom algorithms
using the simpler ``dask.delayed`` interface. This allows one to create graphs
directly with a light annotation of normal python code:

.. code-block:: python

   >>> x = dask.delayed(inc)(1)
   >>> y = dask.delayed(inc)(2)
   >>> z = dask.delayed(add)(x, y)
   >>> z.compute()
   5
   >>> z.visualize()

.. image:: images/inc-add.svg
   :alt: A task graph with two "inc" functions combined using an "add" function resulting in an output node.

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/SHqFmynRxVU"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

Example
-------

Visit https://examples.dask.org/delayed.html to see and run examples using Dask
Delayed.

Sometimes we face problems that are parallelizable, but don't fit into high-level
abstractions like Dask Array or Dask DataFrame.  Consider the following example:

.. code-block:: python

    def inc(x):
        return x + 1

    def double(x):
        return x * 2

    def add(x, y):
        return x + y

    data = [1, 2, 3, 4, 5]

    output = []
    for x in data:
        a = inc(x)
        b = double(x)
        c = add(a, b)
        output.append(c)

    total = sum(output)

There is clearly parallelism in this problem (many of the ``inc``,
``double``, and ``add`` functions can evaluate independently), but it's not
clear how to convert this to a big array or big DataFrame computation.

As written, this code runs sequentially in a single thread.  However, we see that
a lot of this could be executed in parallel.

The Dask ``delayed`` function decorates your functions so that they operate
*lazily*.  Rather than executing your function immediately, it will defer
execution, placing the function and its arguments into a task graph.

.. currentmodule:: dask.delayed

.. autosummary::
    delayed

We slightly modify our code by wrapping functions in ``delayed``.
This delays the execution of the function and generates a Dask graph instead:

.. code-block:: python

    import dask

    output = []
    for x in data:
        a = dask.delayed(inc)(x)
        b = dask.delayed(double)(x)
        c = dask.delayed(add)(a, b)
        output.append(c)

    total = dask.delayed(sum)(output)

We used the ``dask.delayed`` function to wrap the function calls that we want
to turn into tasks.  None of the ``inc``, ``double``, ``add``, or ``sum`` calls
have happened yet. Instead, the object ``total`` is a ``Delayed`` result that
contains a task graph of the entire computation.  Looking at the graph we see
clear opportunities for parallel execution.  The Dask schedulers will exploit
this parallelism, generally improving performance (although not in this
example, because these functions are already very small and fast.)

.. code-block:: python

   total.visualize()  # see image to the right

.. image:: images/delayed-inc-double-add.svg
   :align: right
   :alt: A task graph with many nodes for "inc" and "double" that combine with "add" nodes. The output of the "add" nodes finally aggregate with a "sum" node.

We can now compute this lazy result to execute the graph in parallel:

.. code-block:: python

   >>> total.compute()
   45

Decorator
---------

It is also common to see the delayed function used as a decorator.  Here is a
reproduction of our original problem as a parallel code:

.. code-block:: python

    import dask

    @dask.delayed
    def inc(x):
        return x + 1

    @dask.delayed
    def double(x):
        return x * 2

    @dask.delayed
    def add(x, y):
        return x + y

    data = [1, 2, 3, 4, 5]

    output = []
    for x in data:
        a = inc(x)
        b = double(x)
        c = add(a, b)
        output.append(c)

    total = dask.delayed(sum)(output)


Real time
---------

Sometimes you want to create and destroy work during execution, launch tasks
from other tasks, etc.  For this, see the :doc:`Futures <futures>` interface.


Best Practices
--------------

For a list of common problems and recommendations see :doc:`Delayed Best
Practices <delayed-best-practices>`.


Indirect Dependencies
---------------------

Sometimes you might find yourself wanting to add a dependency to a task that does
not take the result of that dependency as an input. For example when a task depends
on the side-effect of another task. In these cases you can use
``dask.graph_manipulation.bind``.

.. code-block:: python

    import dask
    from dask.graph_manipulation import bind

    DATA = []

    @dask.delayed
    def inc(x):
        return x + 1

    @dask.delayed
    def add_data(x):
        DATA.append(x)

    @dask.delayed
    def sum_data(x):
        return sum(DATA) + x

    a = inc(1)
    b = add_data(a)
    c = inc(3)
    d = add_data(c)
    e = inc(5)
    f = bind(sum_data, [b, d])(e)
    f.compute()

``sum_data`` will operate on DATA only after both the expected items have
been appended to it. ``bind`` can also be used along with direct dependencies
passed through the function arguments.
:orphan:

.. this page is referenced from the topbar which comes from the theme

Community
=========

Dask is used and developed by individuals at a variety of institutions.  It
sits within the broader Python numeric ecosystem commonly referred to as PyData
or SciPy.

Discussion
----------

Conversation happens in the following places:

#.  **Usage questions, requests for help, and general discussions** happen in the
    `Dask Discourse forum`_. If your discussion topic is not a bug report
    or a feature request, this is the best place to start. It's also a good
    place to show off cool things you have built using Dask and to get to know other
    community members.
#.  **Usage questions** may also be directed to `Stack Overflow with the #dask tag`_,
    which is monitored by Dask developers. However, the scope of what is considered
    a good Stack Overflow question can be narrow, so the Dask Discourse forum may
    be a better place to start.
#.  **Bug reports and feature requests** are managed on the `GitHub issue
    tracker`_
#.  **Real-time chat** occurs on
    `https://dask.slack.com/ <https://join.slack.com/t/dask/shared_invite/zt-mfmh7quc-nIrXL6ocgiUH2haLYA914g>`_.
    Note that Slack chat not easily searchable and indexed by search engines, so
    detailed discussion topics around bug reports or usage should go to GitHub issues or
    the Dask Discourse forum, respectively.
#.  **Monthly developer meeting** happens the first Thursday of the month at
    10:00 US Central Time in `this video meeting <https://us06web.zoom.us/j/87619866741?pwd=S2RxMlRKcnVvakt4NHZoS1cwOGZoZz09>`_.
    Meeting notes are available in
    `this Google doc <https://docs.google.com/document/d/1UqNAP87a56ERH_xkQsS5Q_0PKYybd5Lj2WANy_hRzI0/edit>`_.

    .. raw:: html

       <iframe id="calendariframe" src="https://calendar.google.com/calendar/embed?ctz=local&amp;src=4l0vts0c1cgdbq5jhcogj55sfs%40group.calendar.google.com" style="border: 0" width="800" height="600" frameborder="0" scrolling="no"></iframe>
       <script>document.getElementById("calendariframe").src = document.getElementById("calendariframe").src.replace("ctz=local", "ctz=" + Intl.DateTimeFormat().resolvedOptions().timeZone)</script>

    You can subscribe to this calendar to be notified of changes:

    * `Google Calendar <https://calendar.google.com/calendar/u/0?cid=NGwwdnRzMGMxY2dkYnE1amhjb2dqNTVzZnNAZ3JvdXAuY2FsZW5kYXIuZ29vZ2xlLmNvbQ>`__
    * `iCal <https://calendar.google.com/calendar/ical/4l0vts0c1cgdbq5jhcogj55sfs%40group.calendar.google.com/public/basic.ics>`__

.. _`Dask Discourse forum`: https://dask.discourse.group
.. _`Stack Overflow with the #dask tag`: https://stackoverflow.com/questions/tagged/dask
.. _`GitHub issue tracker`: https://github.com/dask/dask/issues/


Asking for help
---------------

We welcome usage questions and bug reports from all users, even those who are
new to using the project.  There are a few things you can do to improve the
likelihood of quickly getting a good answer.

1.  **Ask questions in the right place**:  We strongly prefer the use
    of Discourse or GitHub issues over Slack chat.  Discourse and
    GitHub are more easily searchable by future users, and therefore can be
    useful to many more people than those directly involved.

    If you have a general question about how something should work or
    want best practices then use Discourse.  If you think you have found a
    bug then use GitHub

2.  **Ask only in one place**: Please restrict yourself to posting your
    question in only one place (likely the Dask Discourse or GitHub) and don't post
    in both

3.  **Create a minimal example**:  It is ideal to create `minimal, complete,
    verifiable examples <https://stackoverflow.com/help/mcve>`_.  This
    significantly reduces the time that answerers spend understanding your
    situation, resulting in higher quality answers more quickly.

    See also `this blogpost
    <http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports>`_
    about crafting minimal bug reports.  These have a much higher likelihood of
    being answered


Paid support
------------
In addition to the previous options, paid support is available from the
following organizations (listed in alphabetical order):

-   Anaconda: `<https://www.anaconda.com/products/professional-services>`_
-   Coiled: `<https://coiled.io>`_
-   Quansight: `<https://www.quansight.com/open-source-support>`_
Working with Collections
========================

Often we want to do a bit of custom work with ``dask.delayed`` (for example,
for complex data ingest), then leverage the algorithms in ``dask.array`` or
``dask.dataframe``, and then switch back to custom work.  To this end, all
collections support ``from_delayed`` functions and ``to_delayed``
methods.

As an example, consider the case where we store tabular data in a custom format
not known by Dask DataFrame.  This format is naturally broken apart into
pieces and we have a function that reads one piece into a Pandas DataFrame.
We use ``dask.delayed`` to lazily read these files into Pandas DataFrames,
use ``dd.from_delayed`` to wrap these pieces up into a single
Dask DataFrame, use the complex algorithms within the DataFrame
(groupby, join, etc.), and then switch back to ``dask.delayed`` to save our results
back to the custom format:

.. code-block:: python

   import dask.dataframe as dd
   from dask.delayed import delayed

   from my_custom_library import load, save

   filenames = ...
   dfs = [delayed(load)(fn) for fn in filenames]

   df = dd.from_delayed(dfs)
   df = ... # do work with dask.dataframe

   dfs = df.to_delayed()
   writes = [delayed(save)(df, fn) for df, fn in zip(dfs, filenames)]

   dd.compute(*writes)

Data science is often complex, and ``dask.delayed`` provides a release valve for
users to manage this complexity on their own, and solve the last mile problem
for custom formats and complex situations.
Understanding Performance
=========================

The first step in making computations run quickly is to understand the costs involved.
In Python we often rely on tools like
the `CProfile module <https://docs.python.org/3/library/profile.html>`_,
`%%prun IPython magic <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-prun>`_,
`VMProf <https://vmprof.readthedocs.io/en/latest/>`_, or
`snakeviz <https://jiffyclub.github.io/snakeviz/>`_
to understand the costs associated with our code.
However, few of these tools work well on multi-threaded or multi-process code,
and fewer still on computations distributed among many machines.
We also have new costs like data transfer, serialization, task scheduling overhead, and more
that we may not be accustomed to tracking.

Fortunately, the Dask schedulers come with diagnostics
to help you understand the performance characteristics of your computations.
By using these diagnostics and with some thought,
we can often identify the slow parts of troublesome computations.

The :doc:`single-machine and distributed schedulers <scheduling>` come with *different* diagnostic tools.
These tools are deeply integrated into each scheduler,
so a tool designed for one will not transfer over to the other.

These pages provide four options for profiling parallel code:

1.  :doc:`Visualize task graphs <graphviz>`
2.  :ref:`Single threaded scheduler and a normal Python profiler <single-threaded-scheduler>`
3.  :doc:`Diagnostics for the single-machine scheduler <diagnostics-local>`
4.  :doc:`Diagnostics for the distributed scheduler and dashboard <diagnostics-distributed>`

Additionally, if you are interested in understanding the various phases where
slowdown can occur, you may wish to read the following:

-  :doc:`Phases of computation <phases-of-computation>`
:orphan:

Comparison to Spark
===================

`Apache Spark <https://spark.apache.org/>`_ is a popular distributed computing
tool for tabular datasets that is growing to become a dominant name in Big Data
analysis today.  Dask has several elements that appear to intersect this space
and we are often asked, "How does Dask compare with Spark?"

Answering such comparison questions in an unbiased and informed way is hard,
particularly when the differences can be somewhat technical.  This document
tries to do this; we welcome any corrections.

Summary
-------

Generally Dask is smaller and lighter weight than Spark.  This means that it
has fewer features and, instead, is used in conjunction with other libraries,
particularly those in the numeric Python ecosystem.  It couples with libraries
like Pandas or Scikit-Learn to achieve high-level functionality.

Language
~~~~~~~~

-   Spark is written in Scala with some support for Python and R.  It
    interoperates well with other JVM code.

-   Dask is written in Python and only really supports Python.  It
    interoperates well with C/C++/Fortran/LLVM or other natively compiled
    code linked through Python.

Ecosystem
~~~~~~~~~

-   Spark is an all-in-one project that has inspired its own ecosystem.  It
    integrates well with many other Apache projects.

-   Dask is a component of the larger Python ecosystem.  It couples with and
    enhances other libraries like NumPy, Pandas, and Scikit-Learn.


Age and Trust
~~~~~~~~~~~~~

-   Spark is older (since 2010) and has become a dominant and
    well-trusted tool in the Big Data enterprise world.

-   Dask is younger (since 2014) and is an extension of the
    well trusted NumPy/Pandas/Scikit-learn/Jupyter stack.

Scope
~~~~~

-   Spark is more focused on traditional business intelligence
    operations like SQL and lightweight machine learning.

-   Dask is applied more generally both to business intelligence
    applications, as well as a number of scientific and custom situations.

Internal Design
~~~~~~~~~~~~~~~

-   Spark's internal model is higher level, providing good high level
    optimizations on uniformly applied computations, but lacking flexibility
    for more complex algorithms or ad-hoc systems.  It is fundamentally an
    extension of the Map-Shuffle-Reduce paradigm.

-   Dask's internal model is lower level, and so lacks high level
    optimizations, but is able to implement more sophisticated algorithms and
    build more complex bespoke systems.  It is fundamentally based on generic
    task scheduling.

Scale
~~~~~

-  Spark scales from a single node to thousand-node clusters.
-  Dask scales from a single node to thousand-node clusters.

APIs
~~~~

DataFrames
``````````

-   Spark DataFrame has its own API and memory model.  It also
    implements a large subset of the SQL language.  Spark includes a
    high-level query optimizer for complex queries.

-   Dask DataFrame reuses the Pandas API and memory model.  It implements
    neither SQL nor a query optimizer.  It is able to do random access,
    efficient time series operations, and other Pandas-style indexed
    operations.

Machine Learning
````````````````

-   Spark MLLib is a cohesive project with support for common operations
    that are easy to implement with Spark's Map-Shuffle-Reduce style
    system.  People considering MLLib might also want to consider *other*
    JVM-based machine learning libraries like H2O, which may have better
    performance.

-   Dask relies on and interoperates with existing libraries like
    Scikit-Learn and XGBoost.  These can be more familiar or higher
    performance, but generally results in a less-cohesive whole.  See the
    `dask-ml`_ project for integrations.

Arrays
``````

-   Spark does not include support for multi-dimensional arrays natively
    (this would be challenging given their computation model), although
    some support for two-dimensional matrices may be found in MLLib.
    People may also want to look at the
    `Thunder <https://github.com/thunder-project/thunder>`_ project, which
    combines Apache Spark with NumPy arrays.

-   Dask fully supports the NumPy model for
    :doc:`scalable multi-dimensional arrays <array>`.

Streaming
`````````

-   Spark's support for streaming data is first-class and integrates well
    into their other APIs.  It follows a mini-batch approach.  This
    provides decent performance on large uniform streaming operations.

-   Dask provides a :doc:`real-time futures interface <futures>` that is
    lower-level than Spark streaming.  This enables more creative and
    complex use-cases, but requires more work than Spark streaming.

Graphs / complex networks
`````````````````````````

-  Spark provides GraphX, a library for graph processing.

-  Dask provides no such library.

Custom parallelism
``````````````````

-   Spark generally expects users to compose computations out of their
    high-level primitives (map, reduce, groupby, join, ...).  It is also
    possible to extend Spark through subclassing RDDs, although this is
    rarely done.

-   Dask allows you to specify arbitrary task graphs for more complex and
    custom systems that are not part of the standard set of collections.

.. _dask-ml: https://ml.dask.org


Reasons you might choose Spark
------------------------------

-  You prefer Scala or the SQL language
-  You have mostly JVM infrastructure and legacy systems
-  You want an established and trusted solution for business
-  You are mostly doing business analytics with some lightweight machine learning
-  You want an all-in-one solution


Reasons you might choose Dask
-----------------------------

-  You prefer Python or native code, or have large legacy code bases that you
   do not want to entirely rewrite
-  Your use case is complex or does not cleanly fit the Spark computing model
-  You want a lighter-weight transition from local computing to cluster
   computing
-  You want to interoperate with other technologies and don't mind installing
   multiple packages


Reasons to choose both
----------------------

It is easy to use both Dask and Spark on the same data and on the same cluster.

They can both read and write common formats, like CSV, JSON, ORC, and Parquet,
making it easy to hand results off between Dask and Spark workflows.

They can both deploy on the same clusters.
Most clusters are designed to support many different distributed systems at the
same time, using resource managers like Kubernetes and YARN.  If you already
have a cluster on which you run Spark workloads, it's likely easy to also run
Dask workloads on your current infrastructure and vice versa.

In particular, for users coming from traditional Hadoop/Spark clusters (such as
those sold by Cloudera/Hortonworks) you are using the Yarn resource
manager.  You can deploy Dask on these systems using the `Dask Yarn
<https://yarn.dask.org>`_ project, as well as other projects, like `JupyterHub
on Hadoop <https://jupyterhub-on-hadoop.readthedocs.io/en/latest/>`_.


Developer-Facing Differences
----------------------------

Graph Granularity
~~~~~~~~~~~~~~~~~

Both Spark and Dask represent computations with directed acyclic graphs.  These
graphs however represent computations at very different granularities.

One operation on a Spark RDD might add a node like ``Map`` and ``Filter`` to
the graph.  These are high-level operations that convey meaning and will
eventually be turned into many little tasks to execute on individual workers.
This many-little-tasks state is only available internally to the Spark
scheduler.

Dask graphs skip this high-level representation and go directly to the
many-little-tasks stage.  As such, one ``map`` operation on a Dask collection
will immediately generate and add possibly thousands of tiny tasks to the Dask
graph.

This difference in the scale of the underlying graph has implications on the
kinds of analysis and optimizations one can do and also on the generality that
one exposes to users.  Dask is unable to perform some optimizations that Spark
can because Dask schedulers do not have a top-down picture of the computation
they were asked to perform.  However, Dask is able to easily represent far more
`complex algorithms`_ and expose the creation of these algorithms to normal users.


Conclusion
----------

-   Spark is mature and all-inclusive.  If you want a single project that does
    everything and you're already on Big Data hardware, then Spark is a safe bet,
    especially if your use cases are typical ETL + SQL and you're already using
    Scala.

-   Dask is lighter weight and is easier to integrate into existing code and hardware.
    If your problems vary beyond typical ETL + SQL and you want to add flexible
    parallelism to existing solutions, then Dask may be a good fit, especially if
    you are already using Python and associated libraries like NumPy and Pandas.

If you are looking to manage a terabyte or less of tabular CSV or JSON data,
then you should forget both Spark and Dask and use Postgres_ or MongoDB_.


.. _Spark: https://spark.apache.org/
.. _PySpark: https://spark.apache.org/docs/latest/api/python/
.. _Postgres: https://www.postgresql.org/
.. _MongoDB: https://www.mongodb.org/
.. _`complex algorithms`: http://matthewrocklin.com/blog/work/2015/06/26/Complex-Graphs
Dask Dataframe and SQL
======================

SQL is a method for executing tabular computation on database servers.
Similar operations can be done on Dask Dataframes. Users commonly wish
to link the two together.

This document describes the connection between Dask and SQL-databases
and serves to clarify several of the questions that we commonly
receive from users.

.. contents::
    :local:
    :depth: 1
    :backlinks: top

Does Dask implement SQL?
------------------------

The short answer is "no". Dask has no parser or query planner for SQL
queries. However, the Pandas API, which is largely identical for
Dask Dataframes, has many analogues to SQL operations. A good
description for mapping SQL onto Pandas syntax can be found in the
`pandas docs`_.

.. _pandas docs: https://pandas.pydata.org/docs/getting_started/comparison/comparison_with_sql.html

The following packages may be of interest

- `blazingSQL`_, part of the Rapids project, implements SQL queries using ``cuDF``
  and Dask, for execution on CUDA/GPU-enabled hardware, including referencing
  externally-stored data.

- `dask-sql`_ adds a SQL query layer on top of Dask.
  The API matches blazingSQL but it uses CPU instead of GPU. It still under development
  and not ready for a production use-case.

- `fugue-sql`_ adds an abstract layer that makes code portable between across differing
  computing frameworks such as Pandas, Spark and Dask.

- `pandasql`_ allows executing SQL queries on a pandas table by writing the data to
  ``SQLite``, which may be useful for small toy examples (this package has not been
  maintained for some time).

.. _blazingSQL: https://docs.blazingsql.com/
.. _dask-sql: https://dask-sql.readthedocs.io/en/latest/
.. _fugue-sql: https://fugue-tutorials.readthedocs.io/en/latest/tutorials/fugue_sql/index.html
.. _pandasql: https://github.com/yhat/pandasql/

Database or Dask?
-----------------

A database server is able to process tabular data and produce results just like
Dask Dataframe. Why would you choose to use one over the other?

These days a database server can be a sharded/distributed system, capable of
handling tables with millions of rows. Most database implementations are
geared towards row-wise retrieval and (atomic) updates of small subset of a
table. Configuring a database to be fast for a particular
sort of query can be challenging, but assuming all your data is already in the
database, it may well be the best solution - particularly if you understand
something about SQL query plan optimisation. A SQL implementation can
very efficiently analyse a query to only extract a small part of a table
for consideration, when the rest is excluded by conditionals.

Dask is much more flexible than a database, and designed explicitly
to work with larger-than-memory datasets, in parallel, and potentially distributed
across a cluster. If your workflow is not well suited to SQL, use dask. If
your database server struggles with volume, dask may do better. It
would be best to profile your queries
(and keep in mind other users of the resources!). If you need
to combine data from different sources, dask may be your best option.

You may find the dask API easier to use than writing SQL (if you
are already used to Pandas), and the diagnostic feedback more useful.
These points can debatably be in Dask's favour.

Loading from SQL with read_sql_table or read_sql_query
------------------------------------------------------

Dask allows you to build dataframes from SQL tables and queries using the
function :func:`dask.dataframe.read_sql_table` and :func:`dask.dataframe.read_sql_query`,
based on the `Pandas version`_, sharing most arguments, and using SQLAlchemy
for the actual handling of the queries. You may need to install additional
driver packages for your chosen database server.

.. _Pandas version: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_sql_table.html

Since Dask is designed to work with larger-than-memory datasets, or be distributed
on a cluster, the following are the main differences versus Pandas to watch out for

- Dask does not support arbitrary text queries, only whole tables and SQLAlchemy
  `sql expressions`_

- the con argument must be a `URI string`_, not an SQLAlchemy engine/connection

- partitioning information is *required*, which can be as simple as providing
  an index column argument, or can be more explicit (see below)

- the chunksize argument is not used, since the partitioning must be via an
  index column

.. _URI string: https://docs.sqlalchemy.org/en/13/core/engines.html#database-urls
.. _sql expressions: https://docs.sqlalchemy.org/en/13/core/tutorial.html

If you need something more flexible than this, or the
method fails for you (e.g., on type inference), then skip to the next section.

Why the differences
^^^^^^^^^^^^^^^^^^^

Dask is intended to make processing large volumes of data possible, including
potentially distributing that processing across a cluster. For the retrieval of
data from SQL servers, this means that the query must be partitionable: that
each partition can be fetched independently of others and not dependent on
some global state, and that the definitions of the tasks must be serialisable,
i.e., can be represented as a stream of bytes communicated to workers.

The constraints mean that we cannot directly accept SQLAlchemy engines
or connection objects, since they have internal state (buffers, etc.)
that cannot be serialised. A `URI string`_  must be used, which can be
recreated into a fresh engine on the workers.
Similarly, we cannot accommodate chunked queries
which rely on the internal state of a database cursor; nor LIMIT/OFFSET
queries, which are not guaranteed to be repeatable, and involve scanning
the whole query on th server (which is very inefficient).

**If** your data is small enough not to require Dask's out-of-core and/or
distributed capabilities, then you are probably better to use Pandas or SQLAlchemy
directly.

Index Column
^^^^^^^^^^^^

We need a way to turn a single main query into sub-queries for each
partition. For most reasonable database tables, there should be an obvious
column which can be used for partitioning - it is probably numeric,
and should certainly be indexed in the database. The latter condition
is important, since many simultaneous queries will hit your server once
Dask starts to compute.

By providing just a column name for the index argument, you imply that the
column is numeric, and Dask guesses a reasonable partitioning by evenly
splitting the space between minimum and maximum values into ``npartitions``
intervals. You can also provide the max/min that you would like to
consider so that Dask doesn't need to query for these. Alternatively,
you can have Dask fetch the first few row (5 by default) and use
them to guess the typical bytes/row, and base the partitioning size on
this. Needless to say, the results will vary a lot for tables that are
not uncommonly homogenous.

Specific partitioning
^^^^^^^^^^^^^^^^^^^^^

In some cases, you may have a very good idea of how to partition the data,
for example based on a column that has a finite number of unique values
or categories. This enables using string columns, or anything with a
natural ordering, for the index column, not only numerical types.

In this case, you would provide a specific set of ``divisions``,
the start/end values of the index column for each partition. For example,
if a column happened to contain a random ID in hex string format, then you
could specify 16 partitions with

.. code-block:: python

    df = read_sql_table("mytable", divisions=list("0123456789abcdefh"),
                        index_col="hexID")

so the first partition would have IDs with values ``"0" <= hexID < "1"``, i.e.,
leading character "0".

SQLAlchemy expressions
^^^^^^^^^^^^^^^^^^^^^^

Since we only send the database connection URI and not the engine object,
we cannot rely on SQLAlchemy's table class inference and ORM to conduct queries. However, we can
use the "select" `sql expressions`_, which only get formatted into a text query at
the point of execution.

.. code-block:: python

    from sqlalchemy import sql
    number = sql.column("number")
    name = sql.column("name")
    s1 = sql.select([
            number, name, sql.func.length(name).label("lenname")
        ]
        ).select_from(sql.table("test"))
    data = read_sql_query(
        s1, db, npartitions=2, index_col=number
    )

Here we have also demonstrated the use of the function ``length`` to
perform an operation server-side. Note that it is necessary to *label* such
operations, but you can use them for the index column,
so long as it is also
in the set of selected columns. If using for the index/partitioning, the
column should still be indexed in the database, for performance.
One of the most important functions to consider is ``cast`` to specify the
output data type or conversion in the database, if pandas is having
trouble inferring the data type.

You should be warned, that SQLAlchemy expressions take some time to get
used to, and you can practice with Pandas first, reading only the first small
chunk of a query, until things look right. You can find a more complete
object-oriented example in `this gist`_

.. _this gist: https://gist.github.com/quasiben/08a7f291039db2b04c2e28e1a6c21e3b

Load from SQL, manual approaches
--------------------------------

If ``read_sql_table`` is not sufficient for your needs, you can try one of
the following methods.

Delayed functions
^^^^^^^^^^^^^^^^^

Often you know more about your data and server than the generic approach above
allows. Indeed, some database-like servers may simply not be supported by
SQLAlchemy, or provide an alternate API which is better optimised
(`snowflake example`_).

.. _snowflake example: https://www.saturncloud.io/s/snowflake-and-dask/

If you already have a way to fetch data from the database in partitions,
then you can wrap this function in :func:`dask.delayed` and construct a
dataframe this way. It might look something like

.. code-block:: python

   from dask import delayed
   import dask.dataframe as dd

   @delayed
   def fetch_partition(part):
       conn = establish_connection()
       df = fetch_query(base_query.format(part))
       return df.astype(known_types)

    ddf = dd.from_delayed([fetch_partition(part) for part in parts],
                          meta=known_types,
                          divisions=div_from_parts(parts))

Where you must provide your own functions for setting up a connection to the server,
your own query, and a way to format that query to be specific to each partition.
For example, you might have ranges or specific unique values with a WHERE
clause. The ``known_types`` here is used to transform the dataframe partition and provide
a ``meta``, to help for consistency and avoid Dask having to analyse one partition
up front to guess the columns/types; you may also want to explicitly set the index.
Furthermore, it is a good idea to provide
``divisions`` (the start/end of each partition in the index column), if possible,
since you likely know these from the subqueries you are constructing.

Stream via client
^^^^^^^^^^^^^^^^^

In some cases, the workers may not have access to data, but the client does;
or the initial loading time of the data is not important, so long as the
dataset is then held in cluster memory and available for dask-dataframe
queries. It is possible to construct the dataframe by uploading chunks of
data from the client:

See a complete example of how to do this `here`_

.. _here: https://stackoverflow.com/questions/62818473/why-dasks-read-sql-table-requires-a-index-col-parameter/62821858#62821858


Access data files directly
^^^^^^^^^^^^^^^^^^^^^^^^^^

Some database systems such as Apache Hive store their data in a location
and format that may be directly accessible to Dask, such as parquet files
on S3 or HDFS. In cases where your SQL query would read whole datasets and pass
them to Dask, the streaming of data from the database is very likely the
bottleneck, and it's probably faster to read the source data files directly.

Query pushdown?
---------------

If you define a query based on a database table, then only use some columns
of the output, you may expect that Dask is able to tell the database server
to only send some of the table's data. Dask is not currently able to
do this "pushdown" optimisation, and you would need to change your query using
the SQL expression syntax.
We may be able to resolve this in the future (:issue:`6388`).

If the divisions on your dataframe are well defined, then selections on the
index may successfully avoid reading irrelevant partitions.
Categoricals
============

Dask DataFrame divides `categorical data`_ into two types:

- Known categoricals have the ``categories`` known statically (on the ``_meta``
  attribute).  Each partition **must** have the same categories as found on the
  ``_meta`` attribute
- Unknown categoricals don't know the categories statically, and may have
  different categories in each partition.  Internally, unknown categoricals are
  indicated by the presence of ``dd.utils.UNKNOWN_CATEGORIES`` in the
  categories on the ``_meta`` attribute.  Since most DataFrame operations
  propagate the categories, the known/unknown status should propagate through
  operations (similar to how ``NaN`` propagates)

For metadata specified as a description (option 2 above), unknown categoricals
are created.

Certain operations are only available for known categoricals.  For example,
``df.col.cat.categories`` would only work if ``df.col`` has known categories,
since the categorical mapping is only known statically on the metadata of known
categoricals.

The known/unknown status for a categorical column can be found using the
``known`` property on the categorical accessor:

.. code-block:: python

    >>> ddf.col.cat.known
    False

Additionally, an unknown categorical can be converted to known using
``.cat.as_known()``.  If you have multiple categorical columns in a DataFrame,
you may instead want to use ``df.categorize(columns=...)``, which will convert
all specified columns to known categoricals.  Since getting the categories
requires a full scan of the data, using ``df.categorize()`` is more efficient
than calling ``.cat.as_known()`` for each column (which would result in
multiple scans):

.. code-block:: python

    >>> col_known = ddf.col.cat.as_known()  # use for single column
    >>> col_known.cat.known
    True
    >>> ddf_known = ddf.categorize()        # use for multiple columns
    >>> ddf_known.col.cat.known
    True

To convert a known categorical to an unknown categorical, there is also the
``.cat.as_unknown()`` method. This requires no computation as it's just a
change in the metadata.

Non-categorical columns can be converted to categoricals in a few different
ways:

.. code-block:: python

    # astype operates lazily, and results in unknown categoricals
    ddf = ddf.astype({'mycol': 'category', ...})
    # or
    ddf['mycol'] = ddf.mycol.astype('category')

    # categorize requires computation, and results in known categoricals
    ddf = ddf.categorize(columns=['mycol', ...])

Additionally, with Pandas 0.19.2 and up, ``dd.read_csv`` and ``dd.read_table``
can read data directly into unknown categorical columns by specifying a column
dtype as ``'category'``:

.. code-block:: python

    >>> ddf = dd.read_csv(..., dtype={col_name: 'category'})

.. _`categorical data`: https://pandas.pydata.org/pandas-docs/stable/categorical.html

Moreover, with Pandas 0.21.0 and up, ``dd.read_csv`` and ``dd.read_table`` can read
data directly into *known* categoricals by specifying instances of
``pd.api.types.CategoricalDtype``:

.. code-block:: python

    >>> dtype = {'col': pd.api.types.CategoricalDtype(['a', 'b', 'c'])}
    >>> ddf = dd.read_csv(..., dtype=dtype)

If you write and read to parquet, Dask will forget known categories.
This happens because, due to performance concerns, all the categories are
saved in every partition rather than in the parquet metadata.
It is possible to manually load the categories:

.. code-block:: python

    >>> import dask.dataframe as dd
    >>> import pandas as pd
    >>> df = pd.DataFrame(data=list('abcaabbcc'), columns=['col'])
    >>> df.col = df.col.astype('category')
    >>> ddf = dd.from_pandas(df, npartitions=1)
    >>> ddf.col.cat.known
    True
    >>> ddf.to_parquet('tmp')
    >>> ddf2 = dd.read_parquet('tmp')
    >>> ddf2.col.cat.known
    False
    >>> ddf2.col = ddf2.col.cat.set_categories(ddf2.col.head(1).cat.categories)
    >>> ddf2.col.cat.known
    True
Deploying Clusters
==================

.. toctree::
   :maxdepth: 1
   :hidden:

   deploying-python.rst
   deploying-cli.rst
   deploying-ssh.rst
   deploying-docker.rst
   deploying-hpc.rst
   deploying-kubernetes.rst
   deploying-cloud.rst
   deploying-python-advanced.rst

The ``dask.distributed`` scheduler works well on a single machine and scales to many machines
in a cluster. We recommend using ``dask.distributed`` clusters at all scales for the following
reasons:

1.  It provides access to asynchronous API, notably :doc:`Futures <../../futures>`
2.  It provides a diagnostic dashboard that can provide valuable insight on
    performance and progress
3.  It handles data locality with sophistication, and so can be more
    efficient than the multiprocessing scheduler on workloads that require
    multiple processes

This page describes various ways to set up Dask clusters on different hardware, either
locally on your own machine or on a distributed cluster.  If you are just
getting started then you can save this page for later as Dask runs perfectly well on a single machine
without a distributed scheduler. But once you start using Dask in anger you'll find a lot of benefit
both in terms of scaling and debugging by using the distirbuted scheduler.

You can continue reading or watch the screencast below:

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/TQM9zIBzNBo"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

If you import Dask, set up a computation, and call ``compute``, then you
will use the single-machine scheduler by default.  To use the ``dask.distributed``
scheduler you must set up a ``Client``.

.. code-block:: python

   import dask.dataframe as dd
   df = dd.read_csv(...)
   df.x.sum().compute()  # This uses the single-machine scheduler by default

.. code-block:: python

   from dask.distributed import Client
   client = Client(...)  # Connect to distributed cluster and override default
   df.x.sum().compute()  # This now runs on the distributed system

There are many ways to start the distributed scheduler and worker components that your client
needs to connect to. You can run them manually using :doc:`command line tools <deploying-cli>`
but often the most straight forward way is to use a *cluster manager* utility class.

.. code-block:: python

   from dask.distributed import Client, LocalCluster
   cluster = LocalCluster()  # Launches a scheduler and workers locally
   client = Client(cluster)  # Connect to distributed cluster and override default
   df.x.sum().compute()  # This now runs on the distributed system

There are a number of different *cluster managers* available, so you can use
Dask distributed with a range of platforms. These *cluster managers* deploy a scheduler
and the necessary workers as determined by communicating with the *resource manager*.
All *cluster managers* follow the same interface but have platform specific configuration
options. This makes it convenient to switch from your local machine to a remote multi-node
cluster without sacrificing the flexibility of the platform you are deploying on.

`Dask Jobqueue <https://github.com/dask/dask-jobqueue>`_, for example, is a set of
*cluster managers* for HPC users and works with job queueing systems
(in this case, the *resource manager*) such as `PBS <https://en.wikipedia.org/wiki/Portable_Batch_System>`_,
`Slurm <https://en.wikipedia.org/wiki/Slurm_Workload_Manager>`_,
and `SGE <https://en.wikipedia.org/wiki/Oracle_Grid_Engine>`_.
Those workers are then allocated physical hardware resources.

.. code-block:: python

   from dask.distributed import Client
   from dask_jobqueue import PBSCluster
   cluster = PBSCluster()  # Launches a scheduler and workers on HPC via PBS
   client = Client(cluster)  # Connect to distributed cluster and override default
   df.x.sum().compute()  # This now runs on the distributed system

.. figure:: images/dask-cluster-manager.svg
   :scale: 50%

   An overview of cluster management with Dask distributed.

To summarize, you can use the default, single-machine scheduler to use Dask
on your local machine. If you'd like use a cluster *or* simply take advantage
of the :doc:`extensive diagnostics <../diagnostics-distributed>`,
you can use Dask distributed. The following resources explain
in more detail how to set up Dask on a variety of local and distributed hardware:

- Single Machine:
    - :doc:`Default Scheduler <scheduling>`: The no-setup default.
      Uses local threads or processes for larger-than-memory processing
    - :doc:`dask.distributed <deploying-python>`: The sophistication of
      the newer system on a single machine.  This provides more advanced
      features while still requiring almost no setup.
- Distributed computing:
    - `Beginner's Guide to Configuring a Dask distributed Cluster <https://blog.dask.org/2020/07/30/beginners-config>`_
    - `Overview of cluster management options <https://blog.dask.org/2020/07/23/current-state-of-distributed-dask-clusters>`_
    - :doc:`Manual Setup <deploying-cli>`: The command line interface to set up
      ``dask-scheduler`` and ``dask-worker`` processes.  Useful for IT or
      anyone building a deployment solution.
    - :doc:`SSH <deploying-ssh>`: Use SSH to set up Dask across an un-managed
      cluster.
    - :doc:`High Performance Computers <deploying-hpc>`: How to run Dask on
      traditional HPC environments using tools like MPI, or job schedulers like
      SLURM, SGE, TORQUE, LSF, and so on.
    - :doc:`Kubernetes <deploying-kubernetes>`: Deploy Dask with the
      popular Kubernetes resource manager using either Helm or a native deployment.
    - `YARN / Hadoop <https://yarn.dask.org/en/latest/>`_: Deploy
      Dask on YARN clusters, such as are found in traditional Hadoop
      installations.
    - `Dask Gateway <https://gateway.dask.org/>`_ provides a secure,
      multi-tenant server for managing Dask clusters and allows users to launch
      and use Dask clusters in a shared cluster environment.
    - :doc:`Python API (advanced) <deploying-python-advanced>`: Create
      ``Scheduler`` and ``Worker`` objects from Python as part of a distributed
      Tornado TCP application.  This page is useful for those building custom
      frameworks.
    - :doc:`Docker <deploying-docker>` images are available and may be useful
      in some of the solutions above.
    - :doc:`Cloud <deploying-cloud>` for current recommendations on how to
      deploy Dask and Jupyter on common cloud providers like Amazon, Google, or
      Microsoft Azure.
- Hosted / managed Dask clusters (listed in alphabetical order):
    - `Coiled <https://coiled.io/>`_ handles the creation and management of
      Dask clusters on cloud computing environments (AWS, Azure, and GCP).
    - `Saturn Cloud <https://saturncloud.io/>`_ lets users create
      Dask clusters in a hosted platform or within their own AWS accounts.
Custom Graphs
=============

There may be times when you want to do parallel computing but your application
doesn't fit neatly into something like Dask Array or Dask Bag.  In these
cases, you can interact directly with the Dask schedulers.  These schedulers
operate well as standalone modules.

This separation provides a release valve for complex situations and allows
advanced projects to have additional opportunities for parallel execution, even if
those projects have an internal representation for their computations.  As Dask
schedulers improve or expand to distributed memory, code written to use Dask
schedulers will advance as well.

.. _custom-graph-example:

Example
-------

.. figure:: images/pipeline.svg
   :alt: "Dask graph for data pipeline"
   :align: right

As discussed in the :doc:`motivation <graphs>` and :doc:`specification <spec>`
sections, the schedulers take a task graph (which is a dict of tuples of
functions) and a list of desired keys from that graph.

Here is a mocked out example building a graph for a traditional clean and
analyze pipeline:

.. code-block:: python

   def load(filename):
       ...

   def clean(data):
       ...

   def analyze(sequence_of_data):
       ...

   def store(result):
       with open(..., 'w') as f:
           f.write(result)

   dsk = {'load-1': (load, 'myfile.a.data'),
          'load-2': (load, 'myfile.b.data'),
          'load-3': (load, 'myfile.c.data'),
          'clean-1': (clean, 'load-1'),
          'clean-2': (clean, 'load-2'),
          'clean-3': (clean, 'load-3'),
          'analyze': (analyze, ['clean-%d' % i for i in [1, 2, 3]]),
          'store': (store, 'analyze')}

   from dask.multiprocessing import get
   get(dsk, 'store')  # executes in parallel


Related Projects
----------------

The following excellent projects also provide parallel execution:

*  Joblib_
*  Multiprocessing_
*  `IPython Parallel`_
*  `Concurrent.futures`_
*  `Luigi`_

Each library lets you dictate how your tasks relate to each other with various
levels of sophistication.  Each library executes those tasks with some internal
logic.

Dask schedulers differ in the following ways:

1.  You specify the entire graph as a Python dict rather than using a
    specialized API.
2.  You get a variety of schedulers, ranging from a single-machine, single-core
    scheduler to threaded, multi-process, and distributed options.
3.  You benefit from logic to execute the graph in a way that minimizes memory
    footprint with the Dask single-machine schedulers.

But the other projects offer different advantages and different programming
paradigms.  One should inspect all such projects before selecting one.

.. _Joblib: https://joblib.readthedocs.io/en/latest/
.. _Multiprocessing: https://docs.python.org/3/library/multiprocessing.html
.. _`IPython Parallel`: https://ipyparallel.readthedocs.io/en/latest/
.. _`Concurrent.futures`: https://docs.python.org/3/library/concurrent.futures.html
.. _Luigi: https://luigi.readthedocs.io
Command Line
============

This is the most fundamental way to deploy Dask on multiple machines.  In
production environments this process is often automated by some other resource
manager. Hence, it is rare that people need to follow these instructions
explicitly.  Instead, these instructions are useful to help understand what
*cluster managers* and other automated tooling is doing under the hood and to
help users deploy onto platforms that have no automated tools today.

A ``dask.distributed`` network consists of one ``dask-scheduler`` process and
several ``dask-worker`` processes that connect to that scheduler.  These are
normal Python processes that can be executed from the command line.  We launch
the ``dask-scheduler`` executable in one process and the ``dask-worker``
executable in several processes, possibly on different machines.

To accomplish this, launch ``dask-scheduler`` on one node::

   $ dask-scheduler
   Scheduler at:   tcp://192.0.0.100:8786

Then, launch ``dask-worker`` on the rest of the nodes, providing the address to
the node that hosts ``dask-scheduler``::

   $ dask-worker tcp://192.0.0.100:8786
   Start worker at:  tcp://192.0.0.1:12345
   Registered to:    tcp://192.0.0.100:8786

   $ dask-worker tcp://192.0.0.100:8786
   Start worker at:  tcp://192.0.0.2:40483
   Registered to:    tcp://192.0.0.100:8786

   $ dask-worker tcp://192.0.0.100:8786
   Start worker at:  tcp://192.0.0.3:27372
   Registered to:    tcp://192.0.0.100:8786

The workers connect to the scheduler, which then sets up a long-running network
connection back to the worker.  The workers will learn the location of other
workers from the scheduler.


Handling Ports
--------------

The scheduler and workers both need to accept TCP connections on an open port.
By default, the scheduler binds to port ``8786`` and the worker binds to a
random open port.  If you are behind a firewall then you may have to open
particular ports or tell Dask to listen on particular ports with the ``--port``
and ``--worker-port`` keywords.::

   dask-scheduler --port 8000
   dask-worker --dashboard-address 8000 --nanny-port 8001


Nanny Processes
---------------

Dask workers are run within a nanny process that monitors the worker process
and restarts it if necessary.


Diagnostic Web Servers
----------------------

Additionally, Dask schedulers and workers host interactive diagnostic web
servers using `Bokeh <https://docs.bokeh.org>`_.  These are optional, but
generally useful to users.  The diagnostic server on the scheduler is
particularly valuable, and is served on port ``8787`` by default (configurable
with the ``--dashboard-address`` keyword).

For more information about relevant ports, please take a look at the available
:ref:`command line options <worker-scheduler-cli-options>`.

Automated Tools
---------------

There are various mechanisms to deploy these executables on a cluster, ranging
from manually SSH-ing into all of the machines to more automated systems like
SGE/SLURM/Torque or Yarn/Mesos.  Additionally, cluster SSH tools exist to send
the same commands to many machines.  We recommend searching online for "cluster
ssh" or "cssh".


.. _worker-scheduler-cli-options:

CLI Options
-----------

.. note::

   The command line documentation here may differ depending on your installed
   version. We recommend referring to the output of ``dask-scheduler --help``
   and ``dask-worker --help``.

.. click:: distributed.cli.dask_scheduler:main
   :prog: dask-scheduler
   :show-nested:

.. click:: distributed.cli.dask_worker:main
   :prog: dask-worker
   :show-nested:
Docker Images
=============

Example docker images are maintained at https://github.com/dask/dask-docker
and https://hub.docker.com/r/daskdev/ .

Each image installs the full Dask conda environment (including the distributed
scheduler), Numpy, and Pandas on top of a Miniconda installation on top of
a Debian image.

These images are large, around 1GB.

-   ``daskdev/dask``: This a normal debian + miniconda image with the full Dask
    conda package (including the distributed scheduler), Numpy, and Pandas.
    This image is about 1GB in size.

-   ``daskdev/dask-notebook``: This is based on the
    `Jupyter base-notebook image <https://hub.docker.com/r/jupyter/base-notebook/>`_
    and so it is suitable for use both normally as a Jupyter server, and also as
    part of a JupyterHub deployment.  It also includes a matching Dask software
    environment described above.  This image is about 2GB in size.

Example
-------

Here is a simple example on the local host network

.. code-block:: bash

   docker run -it --network host daskdev/dask dask-scheduler  # start scheduler

   docker run -it --network host daskdev/dask dask-worker localhost:8786 # start worker
   docker run -it --network host daskdev/dask dask-worker localhost:8786 # start worker
   docker run -it --network host daskdev/dask dask-worker localhost:8786 # start worker

   docker run -it --network host daskdev/dask-notebook  # start Jupyter server


Extensibility
-------------

Users can mildly customize the software environment by populating the
environment variables ``EXTRA_APT_PACKAGES``, ``EXTRA_CONDA_PACKAGES``, and
``EXTRA_PIP_PACKAGES``.  If these environment variables are set in the container,
they will trigger calls to the following respectively::

   apt-get install $EXTRA_APT_PACKAGES
   conda install $EXTRA_CONDA_PACKAGES
   python -m pip install $EXTRA_PIP_PACKAGES

For example, the following ``conda`` installs the ``joblib`` package into
the Dask worker software environment:

.. code-block:: bash

   docker run -it -e EXTRA_CONDA_PACKAGES="joblib" daskdev/dask dask-worker localhost:8786

Note that using these can significantly delay the container from starting,
especially when using ``apt``, or ``conda`` (``pip`` is relatively fast).

Remember that it is important for software versions to match between Dask
workers and Dask clients.  As a result, it is often useful to include the same
extra packages in both Jupyter and Worker images.

Source
------

Docker files are maintained at https://github.com/dask/dask-docker.
This repository also includes a docker-compose configuration.
Cloud
=====

There are a variety of ways to deploy Dask on cloud providers.
Cloud providers provide managed services,
like VMs, Kubernetes, Yarn, or custom APIs with which Dask can connect easily.
You may want to consider the following options:

1.  A managed Kubernetes service and Dask's
    :doc:`Kubernetes integration <deploying-kubernetes>`.
2.  A managed Yarn service,
    like `Amazon EMR <https://aws.amazon.com/emr/>`_
    or `Google Cloud DataProc <https://cloud.google.com/dataproc/>`_
    and `Dask-Yarn <https://yarn.dask.org>`_.

    Specific documentation for the popular Amazon EMR service can be found
    `here <https://yarn.dask.org/en/latest/aws-emr.html>`_
3.  Directly launching cloud resources such as VMs or containers via a cluster manager with
    `Dask Cloud Provider <https://cloudprovider.dask.org/en/latest/>`_

Cloud Deployment Example
------------------------

Using `Dask Cloud Provider <https://cloudprovider.dask.org/en/latest/>`_ to launch a cluster of
VMs on a platform like `DigitalOcean <https://www.digitalocean.com/>`_ can be as convenient as
launching a local cluster.

.. code-block:: python

    >>> import dask.config

    >>> dask.config.set({"cloudprovider.digitalocean.token": "yourAPItoken"})

    >>> from dask_cloudprovider.digitalocean import DropletCluster

    >>> cluster = DropletCluster(n_workers=1)
    Creating scheduler instance
    Created droplet dask-38b817c1-scheduler
    Waiting for scheduler to run
    Scheduler is running
    Creating worker instance
    Created droplet dask-38b817c1-worker-dc95260d

Many of the cluster managers in Dask Cloud Provider work by launching VMs with a startup script
that pulls down the :doc:`Dask Docker image <deploying-docker>` and runs Dask components within that container.
As with all cluster managers the VM resources, Docker image, etc are all configurable.

You can then connect a client and work with the cluster as if it were on your local machine.

.. code-block:: python

    >>> from dask.distributed import Client

    >>> client = Client(cluster)

Data Access
-----------

You may want to install additional libraries in your Jupyter and worker images
to access the object stores of each cloud:

-  `s3fs <https://s3fs.readthedocs.io/>`_ for Amazon's S3
-  `gcsfs <https://gcsfs.readthedocs.io/>`_ for Google's GCS
-  `adlfs <https://github.com/dask/adlfs/>`_ for Microsoft's ADL

Historical Libraries
--------------------

Dask previously maintained libraries for deploying Dask on
Amazon's EC2 and Google GKE.
Due to sporadic interest,
and churn both within the Dask library and EC2 itself,
these were not well maintained.
They have since been deprecated in favor of the
:doc:`Kubernetes and Helm <kubernetes-helm>` solution.
Kubernetes and Helm
===================

It is easy to launch a Dask cluster and a Jupyter_ notebook server on cloud
resources using Kubernetes_ and Helm_.

.. _Kubernetes: https://kubernetes.io/
.. _Helm: https://helm.sh/
.. _Jupyter: https://jupyter.org/

This is particularly useful when you want to deploy a fresh Python environment
on Cloud services like Amazon Web Services, Google Compute Engine, or
Microsoft Azure.

If you already have Python environments running in a pre-existing Kubernetes
cluster, then you may prefer the :doc:`Kubernetes native<kubernetes-native>`
documentation, which is a bit lighter weight.

Launch Kubernetes Cluster
-------------------------

This document assumes that you have a Kubernetes cluster and Helm installed.

If this is not the case, then you might consider setting up a Kubernetes cluster
on one of the common cloud providers like Google, Amazon, or
Microsoft.  We recommend the first part of the documentation in the guide
`Zero to JupyterHub <https://zero-to-jupyterhub.readthedocs.io/en/latest/>`_
that focuses on Kubernetes and Helm (you do not need to follow all of these
instructions). In particular, you don't need to install JupyterHub.

- `Creating a Kubernetes Cluster <https://zero-to-jupyterhub.readthedocs.io/en/latest/create-k8s-cluster.html>`_
- `Setting up Helm <https://zero-to-jupyterhub.readthedocs.io/en/latest/setup-helm.html>`_

Alternatively, you may want to experiment with Kubernetes locally using
`Minikube <https://kubernetes.io/docs/getting-started-guides/minikube/>`_.

Which Chart is Right for You?
-----------------------------

Dask maintains a Helm chart repository containing various charts for the Dask community
https://helm.dask.org/ .
You will need to add this to your known channels and update your local charts::

   helm repo add dask https://helm.dask.org/
   helm repo update

We provides two Helm charts. The right one to choose depends on whether you're
deploying Dask for a single user or for many users.


================  =====================================================================
Helm Chart        Use Case
================  =====================================================================
``dask/dask``     Single-user deployment with one notebook server and one Dask Cluster.
``dask/daskhub``  Multi-user deployment with JupyterHub and Dask Gateway.
================  =====================================================================

See :ref:`kubernetes-helm.single` or :ref:`kubernetes-helm.multi` for detailed
instructions on deploying either of these.
As you might suspect, deploying ``dask/daskhub`` is a bit more complicated since
there are more components. If you're just deploying for a single user we'd recommend
using ``dask/dask``.

.. _kubernetes-helm.single:

Helm Install Dask for a Single User
-----------------------------------

Once your Kubernetes cluster is ready, you can deploy dask using the Dask Helm_ chart::

   helm install my-dask dask/dask

This deploys a ``dask-scheduler``, several ``dask-worker`` processes, and
also an optional Jupyter server.


Verify Deployment
^^^^^^^^^^^^^^^^^

This might take a minute to deploy.  You can check its status with
``kubectl``::

   kubectl get pods
   kubectl get services

   $ kubectl get pods
   NAME                                  READY     STATUS              RESTARTS    AGE
   bald-eel-jupyter-924045334-twtxd      0/1       ContainerCreating   0            1m
   bald-eel-scheduler-3074430035-cn1dt   1/1       Running             0            1m
   bald-eel-worker-3032746726-202jt      1/1       Running             0            1m
   bald-eel-worker-3032746726-b8nqq      1/1       Running             0            1m
   bald-eel-worker-3032746726-d0chx      0/1       ContainerCreating   0            1m

   $ kubectl get services
   NAME                 TYPE           CLUSTER-IP      EXTERNAL-IP      PORT(S)                      AGE
   bald-eel-jupyter     LoadBalancer   10.11.247.201   35.226.183.149   80:30173/TCP                  2m
   bald-eel-scheduler   LoadBalancer   10.11.245.241   35.202.201.129   8786:31166/TCP,80:31626/TCP   2m
   kubernetes           ClusterIP      10.11.240.1     <none>           443/TCP
   48m

You can use the addresses under ``EXTERNAL-IP`` to connect to your now-running
Jupyter and Dask systems.

Notice the name ``bald-eel``.  This is the name that Helm has given to your
particular deployment of Dask.  You could, for example, have multiple
Dask-and-Jupyter clusters running at once, and each would be given a different
name.  Note that you will need to use this name to refer to your deployment in the future.
Additionally, you can list all active helm deployments with::

   helm list

   NAME            REVISION        UPDATED                         STATUS      CHART           NAMESPACE
   bald-eel        1               Wed Dec  6 11:19:54 2017        DEPLOYED    dask-0.1.0      default


Connect to Dask and Jupyter
^^^^^^^^^^^^^^^^^^^^^^^^^^^

When we ran ``kubectl get services``, we saw some externally visible IPs:

.. code-block:: bash

   mrocklin@pangeo-181919:~$ kubectl get services
   NAME                 TYPE           CLUSTER-IP      EXTERNAL-IP      PORT(S)                       AGE
   bald-eel-jupyter     LoadBalancer   10.11.247.201   35.226.183.149   80:30173/TCP                  2m
   bald-eel-scheduler   LoadBalancer   10.11.245.241   35.202.201.129   8786:31166/TCP,80:31626/TCP   2m
   kubernetes           ClusterIP      10.11.240.1     <none>           443/TCP                       48m

We can navigate to these services from any web browser. Here, one is the Dask diagnostic
dashboard, and the other is the Jupyter server.  You can log into the Jupyter
notebook server with the password, ``dask``.

You can create a notebook and create a Dask client from there.  The
``DASK_SCHEDULER_ADDRESS`` environment variable has been populated with the
address of the Dask scheduler.  This is available in Python from the ``dask.config`` object.

.. code-block:: python

   >>> import dask
   >>> dask.config.get('scheduler_address')
   'bald-eel-scheduler:8786'

Although you don't need to use this address, the Dask client will find this
variable automatically.

.. code-block:: python

   from dask.distributed import Client, config
   client = Client()


Configure Environment
^^^^^^^^^^^^^^^^^^^^^

By default, the Helm deployment launches three workers using one core each and
a standard conda environment. We can customize this environment by creating a
small yaml file that implements a subset of the values in the
`dask helm chart values.yaml file <https://github.com/dask/helm-chart/blob/main/dask/values.yaml>`_.

For example, we can increase the number of workers, and include extra conda and
pip packages to install on the both the workers and Jupyter server (these two
environments should be matched).

.. code-block:: yaml

   # config.yaml

   worker:
     replicas: 8
     resources:
       limits:
         cpu: 2
         memory: 7.5G
       requests:
         cpu: 2
         memory: 7.5G
     env:
       - name: EXTRA_CONDA_PACKAGES
         value: numba xarray -c conda-forge
       - name: EXTRA_PIP_PACKAGES
         value: s3fs dask-ml --upgrade

   # We want to keep the same packages on the worker and jupyter environments
   jupyter:
     enabled: true
     env:
       - name: EXTRA_CONDA_PACKAGES
         value: numba xarray matplotlib -c conda-forge
       - name: EXTRA_PIP_PACKAGES
         value: s3fs dask-ml --upgrade

This config file overrides the configuration for the number and size of workers and the
conda and pip packages installed on the worker and Jupyter containers.  In
general, we will want to make sure that these two software environments match.

Update your deployment to use this configuration file.  Note that *you will not
use helm install* for this stage: that would create a *new* deployment on the
same Kubernetes cluster.  Instead, you will upgrade your existing deployment by
using the current name::

    helm upgrade bald-eel dask/dask -f config.yaml

This will update those containers that need to be updated.  It may take a minute or so.

As a reminder, you can list the names of deployments you have using ``helm
list``


Check status and logs
^^^^^^^^^^^^^^^^^^^^^

For standard issues, you should be able to see the worker status and logs using the
Dask dashboard (in particular, you can see the worker links from the ``info/`` page).
However, if your workers aren't starting, you can check the status of pods and
their logs with the following commands:

.. code-block:: bash

   kubectl get pods
   kubectl logs <PODNAME>

.. code-block:: bash

   mrocklin@pangeo-181919:~$ kubectl get pods
   NAME                                  READY     STATUS    RESTARTS   AGE
   bald-eel-jupyter-3805078281-n1qk2     1/1       Running   0          18m
   bald-eel-scheduler-3074430035-cn1dt   1/1       Running   0          58m
   bald-eel-worker-1931881914-1q09p      1/1       Running   0          18m
   bald-eel-worker-1931881914-856mm      1/1       Running   0          18m
   bald-eel-worker-1931881914-9lgzb      1/1       Running   0          18m
   bald-eel-worker-1931881914-bdn2c      1/1       Running   0          16m
   bald-eel-worker-1931881914-jq70m      1/1       Running   0          17m
   bald-eel-worker-1931881914-qsgj7      1/1       Running   0          18m
   bald-eel-worker-1931881914-s2phd      1/1       Running   0          17m
   bald-eel-worker-1931881914-srmmg      1/1       Running   0          17m

   mrocklin@pangeo-181919:~$ kubectl logs bald-eel-worker-1931881914-856mm
   EXTRA_CONDA_PACKAGES environment variable found.  Installing.
   Fetching package metadata ...........
   Solving package specifications: .
   Package plan for installation in environment /opt/conda/envs/dask:
   The following NEW packages will be INSTALLED:
       fasteners: 0.14.1-py36_2 conda-forge
       monotonic: 1.3-py36_0    conda-forge
       zarr:      2.1.4-py36_0  conda-forge
   Proceed ([y]/n)?
   monotonic-1.3- 100% |###############################| Time: 0:00:00  11.16 MB/s
   fasteners-0.14 100% |###############################| Time: 0:00:00 576.56 kB/s
   ...


Delete a Helm deployment
^^^^^^^^^^^^^^^^^^^^^^^^

You can always delete a helm deployment using its name::

   helm delete bald-eel --purge

Note that this does not destroy any clusters that you may have allocated on a
Cloud service (you will need to delete those explicitly).


Avoid the Jupyter Server
^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes you do not need to run a Jupyter server alongside your Dask cluster.

.. code-block:: yaml

   jupyter:
     enabled: false

.. _kubernetes-helm.multi:

Helm Install Dask for Mulitple Users
------------------------------------

The ``dask/daskhub`` Helm Chart deploys JupyterHub_, `Dask Gateway`_, and configures
the two to work well together. In particular, Dask Gateway is registered as
a JupyterHub service so that Dask Gateway can re-use JupyterHub's authentication,
and the JupyterHub environment is configured to connect to the Dask Gateway
without any arguments.

.. note::

   The ``dask/daskhub`` helm chart came out of the `Pangeo`_ project, a community
   platform for big data geoscience.

.. _Pangeo: http://pangeo.io/
.. _Dask Gateway: https://gateway.dask.org/
.. _JupyterHub: https://jupyterhub.readthedocs.io/en/stable/

The ``dask/daskhub`` helm chart uses the JupyterHub and Dask-Gateway helm charts.
You'll want to consult the `JupyterHub helm documentation <https://zero-to-jupyterhub.readthedocs.io/en/latest/setup-jupyterhub/setup-jupyterhub.html>`_ and
and `Dask Gateway helm documentation <https://gateway.dask.org/install-kube.html>`_ for further customization. The default values
are at https://github.com/dask/helm-chart/blob/main/daskhub/values.yaml.

Verify that you've set up a Kubernetes cluster and added Dask's helm charts:

.. code-block:: console

   $ helm repo add dask https://helm.dask.org/
   $ helm repo update

JupyterHub and Dask Gateway require a few secret tokens. We'll generate them
on the command line and insert the tokens in a ``secrets.yaml`` file that will
be passed to Helm.

Run the following command, and copy the output. This is our `token-1`.

.. code-block:: console

   $ openssl rand -hex 32  # generate token-1

Run command again and copy the output again. This is our `token-2`.

.. code-block:: console

   $ openssl rand -hex 32  # generate token-2

Now substitute those two values for ``<token-1>`` and ``<token-2>`` below.
Note that ``<token-2>`` is used twice, once for ``jupyterhub.hub.services.dask-gateway.apiToken``, and a second time for ``dask-gateway.gateway.auth.jupyterhub.apiToken``.

.. code-block:: yaml

   # file: secrets.yaml
   jupyterhub:
     proxy:
       secretToken: "<token-1>"
     hub:
       services:
         dask-gateway:
           apiToken: "<token-2>"

   dask-gateway:
     gateway:
       auth:
         jupyterhub:
           apiToken: "<token-2>"

Now we're ready to install DaskHub

.. code-block:: console

   $ helm upgrade --wait --install --render-subchart-notes \
       dhub dask/daskhub \
       --values=secrets.yaml


The output explains how to find the IPs for your JupyterHub depoyment.

.. code-block:: console

   $ kubectl get service proxy-public
   NAME           TYPE           CLUSTER-IP      EXTERNAL-IP      PORT(S)                      AGE
   proxy-public   LoadBalancer   10.43.249.239   35.202.158.223   443:31587/TCP,80:30500/TCP   2m40s


Creating a Dask Cluster
^^^^^^^^^^^^^^^^^^^^^^^

To create a Dask cluster on this deployment, users need to connect to the Dask Gateway

.. code-block:: python

   >>> from dask_gateway import GatewayCluster
   >>> cluster = GatewayCluster()
   >>> client = cluster.get_client()
   >>> cluster

Depending on the configuration, users may need to ``cluster.scale(n)`` to
get workers. See https://gateway.dask.org/ for more on Dask Gateway.

Matching the User Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Dask Clients will be running the JupyterHub's singleuser environment. To ensure
that the same environment is used for the scheduler and workers, you can provide
it as a Gateway option and configure the ``singleuser`` environment to default
to the value set by JupyterHub.

.. code-block:: yaml

   # config.yaml
   jupyterhub:
     singleuser:
       extraEnv:
         DASK_GATEWAY__CLUSTER__OPTIONS__IMAGE: '{JUPYTER_IMAGE_SPEC}'

   dask-gateway:
     gateway:
       extraConfig:
         optionHandler: |
           from dask_gateway_server.options import Options, Integer, Float, String
           def option_handler(options):
               if ":" not in options.image:
                   raise ValueError("When specifying an image you must also provide a tag")
               return {
                   "image": options.image,
               }
           c.Backend.cluster_options = Options(
               String("image", default="pangeo/base-notebook:2020.07.28", label="Image"),
               handler=option_handler,
           )

The user environment will need to include ``dask-gateway``. Any packages installed
manually after the ``singleuser`` pod started will not be included in the worker
environment.
Extending DataFrames
====================

Subclass DataFrames
-------------------

There are a few projects that subclass or replicate the functionality of Pandas
objects:

-  GeoPandas: for Geospatial analytics
-  cuDF: for data analysis on GPUs
-  ...

These projects may also want to produce parallel variants of themselves with
Dask, and may want to reuse some of the code in Dask DataFrame. Subclassing
Dask DataFrames is intended for maintainers of these libraries and not for
general users.


Implement dask, name, meta, and divisions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You will need to implement ``._meta``, ``.dask``, ``.divisions``, and
``._name`` as defined in the :doc:`DataFrame design docs <dataframe-design>`.


Extend Dispatched Methods
^^^^^^^^^^^^^^^^^^^^^^^^^

If you are going to pass around Pandas-like objects that are not normal Pandas
objects, then we ask you to extend a few dispatched methods: ``make_meta``,
``get_parallel_type``, and ``concat``.

make_meta
"""""""""

This function returns an empty version of one of your non-Dask objects, given a
non-empty non-Dask object:

.. code-block:: python

   from dask.dataframe.dispatch import make_meta_dispatch

   @make_meta_dispatch.register(MyDataFrame)
   def make_meta_dataframe(df, index=None):
       return df.head(0)


   @make_meta_dispatch.register(MySeries)
   def make_meta_series(s, index=None):
       return s.head(0)


   @make_meta_dispatch.register(MyIndex)
   def make_meta_index(ind, index=None):
       return ind[:0]

For dispatching any arbitrary ``object`` types to a respective back-end, we
recommend registering a dispatch for ``make_meta_obj``:

.. code-block:: python

    from dask.dataframe.dispatch import make_meta_obj

    @make_meta_obj.register(MyDataFrame)
    def make_meta_object(x, index=None):
        if isinstance(x, dict):
            return MyDataFrame()
        elif isinstance(x, int):
            return MySeries
        .
        .
        .


Additionally, you should create a similar function that returns a non-empty
version of your non-Dask DataFrame objects filled with a few rows of
representative or random data.  This is used to guess types when they are not
provided.  It should expect an empty version of your object with columns,
dtypes, index name, and it should return a non-empty version:

.. code-block:: python

   from dask.dataframe.utils import meta_nonempty

   @meta_nonempty.register(MyDataFrame)
   def meta_nonempty_dataframe(df):
       ...
       return MyDataFrame(..., columns=df.columns,
                          index=MyIndex(..., name=df.index.name)


   @meta_nonempty.register(MySeries)
   def meta_nonempty_series(s):
       ...


   @meta_nonempty.register(MyIndex)
   def meta_nonempty_index(ind):
       ...


get_parallel_type
"""""""""""""""""

Given a non-Dask DataFrame object, return the Dask equivalent:

.. code-block:: python

   from dask.dataframe.core import get_parallel_type

   @get_parallel_type.register(MyDataFrame)
   def get_parallel_type_dataframe(df):
       return MyDaskDataFrame


   @get_parallel_type.register(MySeries)
   def get_parallel_type_series(s):
       return MyDaskSeries


   @get_parallel_type.register(MyIndex)
   def get_parallel_type_index(ind):
       return MyDaskIndex


concat
""""""

Concatenate many of your non-Dask DataFrame objects together.  It should expect
a list of your objects (homogeneously typed):

.. code-block:: python

   from dask.dataframe.methods import concat_dispatch

   @concat_dispatch.register((MyDataFrame, MySeries, MyIndex))
   def concat_pandas(dfs, axis=0, join='outer', uniform=False, filter_warning=True):
       ...


.. _extensionarrays:

Extension Arrays
----------------

Rather than subclassing Pandas DataFrames, you may be interested in extending
Pandas with `Extension Arrays
<https://pandas.pydata.org/pandas-docs/stable/extending.html>`_.

All of the first-party extension arrays (those implemented in pandas itself)
are supported directly by dask.

Developers implementing third-party extension arrays (outside of pandas) will
need to do register their ``ExtensionDtype`` with Dask so that it works
correctly in ``dask.dataframe``.

For example, we'll register the *test-only* ``DecimalDtype`` from pandas
test suite.

.. code-block:: python

   from decimal import Decimal
   from dask.dataframe.extensions import make_array_nonempty, make_scalar
   from pandas.tests.extension.decimal import DecimalArray, DecimalDtype

   @make_array_nonempty.register(DecimalDtype)
   def _(dtype):
       return DecimalArray._from_sequence([Decimal('0'), Decimal('NaN')],
                                          dtype=dtype)


   @make_scalar.register(Decimal)
   def _(x):
      return Decimal('1')


Internally, Dask will use this to create a small dummy Series for tracking
metadata through operations.

.. code-block:: python

   >>> make_array_nonempty(DecimalDtype())
   <DecimalArray>
   [Decimal('0'), Decimal('NaN')]
   Length: 2, dtype: decimal

So you (or your users) can now create and store a dask ``DataFrame`` or
``Series`` with your extension array contained within.

.. code-block:: python

   >>> from decimal import Decimal
   >>> import dask.dataframe as dd
   >>> import pandas as pd
   >>> from pandas.tests.extension.decimal import DecimalArray
   >>> ser = pd.Series(DecimalArray([Decimal('0.0')] * 10))
   >>> dser = dd.from_pandas(ser, 3)
   >>> dser
   Dask Series Structure:
   npartitions=3
   0    decimal
   4        ...
   8        ...
   9        ...
   dtype: decimal
   Dask Name: from_pandas, 3 tasks

Notice the ``decimal`` dtype.

.. _dataframe.accessors:

Accessors
---------

Many extension arrays expose their functionality on Series or DataFrame objects
using accessors. Dask provides decorators to register accessors similar to pandas. See
`the pandas documentation on accessors <http://pandas.pydata.org/pandas-docs/stable/development/extending.html#registering-custom-accessors>`_
for more.

.. currentmodule:: dask.dataframe

.. autofunction:: dask.dataframe.extensions.register_dataframe_accessor
.. autofunction:: dask.dataframe.extensions.register_series_accessor
.. autofunction:: dask.dataframe.extensions.register_index_accessor
10 Minutes to Dask
==================

This is a short overview of what you can do with Dask. It is geared towards new users.
There is much more information contained in the rest of the documentation.

We normally import dask as follows:

.. code-block:: python

   >>> import numpy as np
   >>> import pandas as pd

   >>> import dask.dataframe as dd
   >>> import dask.array as da
   >>> import dask.bag as db

Based on the type of data you are working with, you might not need all of these.

Create a High-Level Collection
------------------------------

You can make a Dask collection from scratch by supplying existing data and optionally
including information about how the chunks should be structured.

.. tabs::

   .. group-tab:: DataFrame

      .. code-block:: python

         >>> index = pd.date_range("2021-09-01", periods=2400, freq="1H")
         ... df = pd.DataFrame({"a": np.arange(2400), "b": list("abcaddbe" * 300)}, index=index)
         ... ddf = dd.from_pandas(df, npartitions=10)
         ... ddf
         Dask DataFrame Structure:
                                 a       b
         npartitions=10
         2021-09-01 00:00:00  int64  object
         2021-09-11 00:00:00    ...     ...
         ...                    ...     ...
         2021-11-30 00:00:00    ...     ...
         2021-12-09 23:00:00    ...     ...
         Dask Name: from_pandas, 10 tasks

      Now we have a DataFrame with 2 columns and 2400 rows composed of 10 partitions where
      each partition has 240 rows. Each partition represents a piece of the data.

      Here are some key properties of an DataFrame:

      .. code-block:: python

         >>> # check the index values covered by each partition
         ... ddf.divisions
         (Timestamp('2021-09-01 00:00:00', freq='H'),
         Timestamp('2021-09-11 00:00:00', freq='H'),
         Timestamp('2021-09-21 00:00:00', freq='H'),
         Timestamp('2021-10-01 00:00:00', freq='H'),
         Timestamp('2021-10-11 00:00:00', freq='H'),
         Timestamp('2021-10-21 00:00:00', freq='H'),
         Timestamp('2021-10-31 00:00:00', freq='H'),
         Timestamp('2021-11-10 00:00:00', freq='H'),
         Timestamp('2021-11-20 00:00:00', freq='H'),
         Timestamp('2021-11-30 00:00:00', freq='H'),
         Timestamp('2021-12-09 23:00:00', freq='H'))

         >>> # access a particular partition
         ... ddf.partitions[1]
         Dask DataFrame Structure:
                           a       b
         npartitions=1
         2021-09-11     int64  object
         2021-09-21       ...     ...
         Dask Name: blocks, 11 tasks

   .. group-tab:: Array

      .. code-block:: python

         >>> data = np.arange(100_000).reshape(200, 500)
         ... a = da.from_array(data, chunks=(100, 100))
         ... a
         dask.array<array, shape=(200, 500), dtype=int64, chunksize=(100, 100), chunktype=numpy.ndarray>

      Now we have a 2D array with the shape (200, 500) composed of 10 chunks where
      each chunk has the shape (100, 100). Each chunk represents a piece of the data.

      Here are some key properties of an Array:

      .. code-block:: python

         >>> # inspect the chunks
         ... a.chunks
         ((100, 100), (100, 100, 100, 100, 100))

         >>> # access a particular block of data
         ... a.blocks[1, 3]
         dask.array<blocks, shape=(100, 100), dtype=int64, chunksize=(100, 100), chunktype=numpy.ndarray>

   .. group-tab:: Bag

      .. code-block:: python

         >>> b = db.from_sequence([1, 2, 3, 4, 5, 6, 2, 1], npartitions=2)
         ... b
         dask.bag<from_sequence, npartitions=2>

      Now we have a sequence with 8 items composed of 2 partitions where each partition
      has 4 items in it. Each partition represents a piece of the data.


Indexing
--------

Indexing Dask collections feels just like slicing numpy arrays or pandas dataframes.

.. tabs::

   .. group-tab:: DataFrame

      .. code-block:: python

         >>> ddf.b
         Dask Series Structure:
         npartitions=10
         2021-09-01 00:00:00    object
         2021-09-11 00:00:00       ...
                                 ...
         2021-11-30 00:00:00       ...
         2021-12-09 23:00:00       ...
         Name: b, dtype: object
         Dask Name: getitem, 20 tasks

         >>> ddf["2021-10-01": "2021-10-09 5:00"]
         Dask DataFrame Structure:
                                          a       b
         npartitions=1
         2021-10-01 00:00:00.000000000  int64  object
         2021-10-09 05:00:59.999999999    ...     ...
         Dask Name: loc, 11 tasks

   .. group-tab:: Array

      .. code-block:: python

         >>> a[:50, 200]
         dask.array<getitem, shape=(50,), dtype=int64, chunksize=(50,), chunktype=numpy.ndarray>

   .. group-tab:: Bag

      A Bag is an unordered collection allowing repeats. So it is like a list, but it doesn’t
      guarantee an ordering among elements. There is no way to index Bags since they are
      not ordered.


Computation
-----------

Dask is lazily evaluated. The result from a computation isn't computed until
you ask for it. Instead, a Dask task graph for the computation is produced.

Anytime you have a Dask object and you want to get the result, call ``compute``:

.. tabs::

   .. group-tab:: DataFrame

      .. code-block:: python

         >>> ddf["2021-10-01": "2021-10-09 5:00"].compute()
                              a  b
         2021-10-01 00:00:00  720  a
         2021-10-01 01:00:00  721  b
         2021-10-01 02:00:00  722  c
         2021-10-01 03:00:00  723  a
         2021-10-01 04:00:00  724  d
         ...                  ... ..
         2021-10-09 01:00:00  913  b
         2021-10-09 02:00:00  914  c
         2021-10-09 03:00:00  915  a
         2021-10-09 04:00:00  916  d
         2021-10-09 05:00:00  917  d

         [198 rows x 2 columns]

   .. group-tab:: Array

      .. code-block:: python

         >>> a[:50, 200].compute()
         array([  200,   700,  1200,  1700,  2200,  2700,  3200,  3700,  4200,
               4700,  5200,  5700,  6200,  6700,  7200,  7700,  8200,  8700,
               9200,  9700, 10200, 10700, 11200, 11700, 12200, 12700, 13200,
               13700, 14200, 14700, 15200, 15700, 16200, 16700, 17200, 17700,
               18200, 18700, 19200, 19700, 20200, 20700, 21200, 21700, 22200,
               22700, 23200, 23700, 24200, 24700])

   .. group-tab:: Bag

      .. code-block:: python

         >>> b.compute()
         [1, 2, 3, 4, 5, 6, 2, 1]


Methods
-------

Dask collections match existing numpy and pandas methods, so they should feel familiar.
Call the method to set up the task graph, and then call ``compute`` to get the result.

.. tabs::

   .. group-tab:: DataFrame

      .. code-block:: python

         >>> ddf.a.mean()
         dd.Scalar<series-..., dtype=float64>

         >>> ddf.a.mean().compute()
         1199.5

         >>> ddf.b.unique()
         Dask Series Structure:
         npartitions=1
            object
               ...
         Name: b, dtype: object
         Dask Name: unique-agg, 33 tasks

         >>> ddf.b.unique().compute()
         0    a
         1    b
         2    c
         3    d
         4    e
         Name: b, dtype: object

      Methods can be chained together just like in pandas

      .. code-block:: python

         >>> result = ddf["2021-10-01": "2021-10-09 5:00"].a.cumsum() - 100
         ... result
         Dask Series Structure:
         npartitions=1
         2021-10-01 00:00:00.000000000    int64
         2021-10-09 05:00:59.999999999      ...
         Name: a, dtype: int64
         Dask Name: sub, 16 tasks

         >>> result.compute()
         2021-10-01 00:00:00       620
         2021-10-01 01:00:00      1341
         2021-10-01 02:00:00      2063
         2021-10-01 03:00:00      2786
         2021-10-01 04:00:00      3510
                                 ...
         2021-10-09 01:00:00    158301
         2021-10-09 02:00:00    159215
         2021-10-09 03:00:00    160130
         2021-10-09 04:00:00    161046
         2021-10-09 05:00:00    161963
         Freq: H, Name: a, Length: 198, dtype: int64

   .. group-tab:: Array

      .. code-block:: python

         >>> a.mean()
         dask.array<mean_agg-aggregate, shape=(), dtype=float64, chunksize=(), chunktype=numpy.ndarray>

         >>> a.mean().compute()
         49999.5

         >>> np.sin(a)
         dask.array<sin, shape=(200, 500), dtype=float64, chunksize=(100, 100), chunktype=numpy.ndarray>

         >>> np.sin(a).compute()
         array([[ 0.        ,  0.84147098,  0.90929743, ...,  0.58781939,
                  0.99834363,  0.49099533],
               [-0.46777181, -0.9964717 , -0.60902011, ..., -0.89796748,
               -0.85547315, -0.02646075],
               [ 0.82687954,  0.9199906 ,  0.16726654, ...,  0.99951642,
                  0.51387502, -0.4442207 ],
               ...,
               [-0.99720859, -0.47596473,  0.48287891, ..., -0.76284376,
                  0.13191447,  0.90539115],
               [ 0.84645538,  0.00929244, -0.83641393, ...,  0.37178568,
               -0.5802765 , -0.99883514],
               [-0.49906936,  0.45953849,  0.99564877, ...,  0.10563876,
                  0.89383946,  0.86024828]])

         >>> a.T
         dask.array<transpose, shape=(500, 200), dtype=int64, chunksize=(100, 100), chunktype=numpy.ndarray>

         >>> a.T.compute()
         array([[    0,   500,  1000, ..., 98500, 99000, 99500],
               [    1,   501,  1001, ..., 98501, 99001, 99501],
               [    2,   502,  1002, ..., 98502, 99002, 99502],
               ...,
               [  497,   997,  1497, ..., 98997, 99497, 99997],
               [  498,   998,  1498, ..., 98998, 99498, 99998],
               [  499,   999,  1499, ..., 98999, 99499, 99999]])

      Methods can be chained together just like in NumPy

      .. code-block:: python

         >>> b = a.max(axis=1)[::-1] + 10
         ... b
         dask.array<add, shape=(200,), dtype=int64, chunksize=(100,), chunktype=numpy.ndarray>

         >>> b[:10].compute()
         array([100009,  99509,  99009,  98509,  98009,  97509,  97009,  96509,
               96009,  95509])

   .. group-tab:: Bag

      Dask Bag implements operations like ``map``, ``filter``, ``fold``, and
      ``groupby`` on collections of generic Python objects.

      .. code-block:: python

         >>> b.filter(lambda x: x % 2)
         dask.bag<filter-lambda, npartitions=2>

         >>> b.filter(lambda x: x % 2).compute()
         [1, 3, 5, 1]

         >>> b.distinct()
         dask.bag<distinct-aggregate, npartitions=1>

         >>> b.distinct().compute()
         [1, 2, 3, 4, 5, 6]

      Methods can be chained together.

      .. code-block:: python

         >>> c = db.zip(b, b.map(lambda x: x * 10))
         ... c
         dask.bag<zip, npartitions=2>

         >>> c.compute()
         [(1, 10), (2, 20), (3, 30), (4, 40), (5, 50), (6, 60), (2, 20), (1, 10)]


Visualize the Task Graph
------------------------

So far we've been setting up computations and calling ``compute``. In addition to
triggering computation, we can inspect the task graph to figure out what's going on.

.. tabs::

   .. group-tab:: DataFrame

      .. code-block:: python

         >>> result.dask
         HighLevelGraph with 7 layers.
         <dask.highlevelgraph.HighLevelGraph object at 0x7f129df7a9d0>
         0. from_pandas-0b850a81e4dfe2d272df4dc718065116
         1. loc-fb7ada1e5ba8f343678fdc54a36e9b3e
         2. getitem-55d10498f88fc709e600e2c6054a0625
         3. series-cumsum-map-131dc242aeba09a82fea94e5442f3da9
         4. series-cumsum-take-last-9ebf1cce482a441d819d8199eac0f721
         5. series-cumsum-d51d7003e20bd5d2f767cd554bdd5299
         6. sub-fed3e4af52ad0bd9c3cc3bf800544f57

         >>> result.visualize()

      .. image:: images/10_minutes_dataframe_graph.png
         :alt: Dask task graph for the Dask dataframe computation. The task graph shows a "loc" and "getitem" operations selecting a small section of the dataframe values, before applying a cumulative sum "cumsum" operation, then finally subtracting a value from the result.

   .. group-tab:: Array

      .. code-block:: python

         >>> b.dask
         HighLevelGraph with 6 layers.
         <dask.highlevelgraph.HighLevelGraph object at 0x7fd33a4aa400>
         0. array-ef3148ecc2e8957c6abe629e08306680
         1. amax-b9b637c165d9bf139f7b93458cd68ec3
         2. amax-partial-aaf8028d4a4785f579b8d03ffc1ec615
         3. amax-aggregate-07b2f92aee59691afaf1680569ee4a63
         4. getitem-f9e225a2fd32b3d2f5681070d2c3d767
         5. add-f54f3a929c7efca76a23d6c42cdbbe84

         >>> b.visualize()

      .. image:: images/10_minutes_array_graph.png
         :alt: Dask task graph for the Dask array computation. The task graph shows many "amax" operations on each chunk of the Dask array, that are then aggregated to find "amax" along the first array axis, then reversing the order of the array values with a "getitem" slicing operation, before an "add" operation to get the final result.

   .. group-tab:: Bag

      .. code-block:: python

         >>> c.dask
         HighLevelGraph with 3 layers.
         <dask.highlevelgraph.HighLevelGraph object at 0x7f96d0814fd0>
         0. from_sequence-cca2a33ba6e12645a0c9bc0fd3fe6c88
         1. lambda-93a7a982c4231fea874e07f71b4bcd7d
         2. zip-474300792cc4f502f1c1f632d50e0272

         >>> c.visualize()

      .. image:: images/10_minutes_bag_graph.png
         :alt: Dask task graph for the Dask bag computation. The task graph shows a "lambda" operation, and then a "zip" operation is applied to the partitions of the Dask bag. There is no communication needed between the bag partitions, this is an embarrassingly parallel computation.

Low-Level Interfaces
--------------------
Often when parallelizing existing code bases or building custom algorithms, you
run into code that is parallelizable, but isn't just a big DataFrame or array.

.. tabs::

   .. group-tab:: Delayed: Lazy

      Dask Delayed let you to wrap individual function calls into a lazily constructed task graph:

      .. code-block:: python

         import dask

         @dask.delayed
         def inc(x):
            return x + 1

         @dask.delayed
         def add(x, y):
            return x + y

         a = inc(1)       # no work has happened yet
         b = inc(2)       # no work has happened yet
         c = add(a, b)    # no work has happened yet

         c = c.compute()  # This triggers all of the above computations

   .. group-tab:: Futures: Immediate

      Unlike the interfaces described so far, Futures are eager. Computation starts as soon
      as the function is submitted.

      .. code-block:: python

         from dask.distributed import Client

         client = Client()

         def inc(x):
            return x + 1

         def add(x, y):
            return x + y

         a = client.submit(inc, 1)     # work starts immediately
         b = client.submit(inc, 2)     # work starts immediately
         c = client.submit(add, a, b)  # work starts immediately

         c = c.result()                # block until work finishes, then gather result

      .. note::

         Futures can only be used with distributed cluster. See the section below for more
         information.


Scheduling
----------

After you have generated a task graph, it is the scheduler's job to execute it.

By default when you call ``compute`` on a Dask object, Dask uses the thread
pool on your computer to run computations in parallel.

If you want more control, use the distributed scheduler instead. Despite having
"distributed" in it's name, the distributed scheduler works well
on both single and multiple machines. Think of it as the "advanced scheduler".

.. tabs::

   .. group-tab:: Local

      This is how you set up a cluster that uses only your own computer.

      .. code-block:: python

         >>> from dask.distributed import Client
         ...
         ... client = Client()
         ... client
         <Client: 'tcp://127.0.0.1:41703' processes=4 threads=12, memory=31.08 GiB>

   .. group-tab:: Remote

      This is how you connect to a cluster that is already running.

      .. code-block:: python

         >>> from dask.distributed import Client
         ...
         ... client = Client("<url-of-scheduler>")
         ... client
         <Client: 'tcp://127.0.0.1:41703' processes=4 threads=12, memory=31.08 GiB>

      There are a variety of ways to set up a remote cluster. Refer to
      :doc:`how to deploy dask clusters <how-to/deploy-dask-clusters>` for more
      information.

Once you create a client, any computation will run on the cluster that it points to.


Diagnostics
-----------

When using a distributed cluster, Dask provides a diagnostics dashboard where you can
see your tasks as they are processed.

.. code-block:: python

   >>> client.dashboard_link
   'http://127.0.0.1:8787/status'

To learn more about those graphs take a look at :doc:`diagnostics-distributed`.
Best Practices
==============

It is easy to get started with Dask's APIs, but using them *well* requires some
experience. This page contains suggestions for best practices, and includes
solutions to common problems.

This document specifically focuses on best practices that are shared among all
of the Dask APIs.  Readers may first want to investigate one of the
API-specific Best Practices documents first.

-  :doc:`Arrays <array-best-practices>`
-  :doc:`DataFrames <dataframe-best-practices>`
-  :doc:`Delayed <delayed-best-practices>`


Start Small
-----------

Parallelism brings extra complexity and overhead.
Sometimes it's necessary for larger problems, but often it's not.
Before adding a parallel computing system like Dask to your workload you may
want to first try some alternatives:

-   **Use better algorithms or data structures**:  NumPy, Pandas, Scikit-Learn
    may have faster functions for what you're trying to do.  It may be worth
    consulting with an expert or reading through their docs again to find a
    better pre-built algorithm.

-   **Better file formats**:  Efficient binary formats that support random
    access can often help you manage larger-than-memory datasets efficiently and
    simply.  See the `Store Data Efficiently`_ section below.

-   **Compiled code**:  Compiling your Python code with Numba or Cython might
    make parallelism unnecessary.  Or you might use the multi-core parallelism
    available within those libraries.

-   **Sampling**:  Even if you have a lot of data, there might not be much
    advantage from using all of it.  By sampling intelligently you might be able
    to derive the same insight from a much more manageable subset.

-   **Profile**:  If you're trying to speed up slow code it's important that
    you first understand why it is slow.  Modest time investments in profiling
    your code can help you to identify what is slowing you down.  This
    information can help you make better decisions about if parallelism is likely
    to help, or if other approaches are likely to be more effective.


Use The Dashboard
-----------------

Dask's dashboard helps you to understand the state of your workers.
This information can help to guide you to efficient solutions.
In parallel and distributed computing there are new costs to be aware of and so
your old intuition may no longer be true.  Working with the dashboard can help
you relearn about what is fast and slow and how to deal with it.

See :doc:`Documentation on Dask's dashboard <diagnostics-distributed>` for more
information.


Avoid Very Large Partitions
---------------------------

Your chunks of data should be small enough so that many of them fit in a
worker's available memory at once.  You often control this when you select
partition size in Dask DataFrame or chunk size in Dask Array.

Dask will likely manipulate as many chunks in parallel on one machine as you
have cores on that machine.  So if you have 1 GB chunks and ten
cores, then Dask is likely to use *at least* 10 GB of memory.  Additionally,
it's common for Dask to have 2-3 times as many chunks available to work on so
that it always has something to work on.

If you have a machine with 100 GB and 10 cores, then you might want to choose
chunks in the 1GB range.  You have space for ten chunks per core which gives
Dask a healthy margin, without having tasks that are too small

Note that you also want to avoid chunk sizes that are too small.  See the next
section for details.


Avoid Very Large Graphs
-----------------------

Dask workloads are composed of *tasks*.
A task is a Python function, like ``np.sum`` applied onto a Python object,
like a Pandas dataframe or NumPy array.  If you are working with Dask
collections with many partitions, then every operation you do, like ``x + 1``
likely generates many tasks, at least as many as partitions in your collection.

Every task comes with some overhead.  This is somewhere between 200us and 1ms.
If you have a computation with thousands of tasks this is fine, there will be
about a second of overhead, and that may not trouble you.

However when you have very large graphs with millions of tasks then this may
become troublesome, both because overhead is now in the 10 minutes to hours
range, and also because the overhead of dealing with such a large graph can
start to overwhelm the scheduler.

There are a few things you can do to address this:

-   Build smaller graphs.  You can do this by ...

    -  **Increasing your chunk size:**  If you have a 1000 GB of data and are using
       10 MB chunks, then you have 100,000 partitions.  Every operation on such
       a collection will generate at least 100,000 tasks.

       However if you increase your chunksize to 1 GB or even a few GB then you
       reduce the overhead by orders of magnitude.  This requires that your
       workers have much more than 1 GB of memory, but that's typical for larger
       workloads.

    -  **Fusing operations together:** Dask will do a bit of this on its own, but you
       can help it.  If you have a very complex operation with dozens of
       sub-operations, maybe you can pack that into a single Python function
       and use a function like ``da.map_blocks`` or ``dd.map_partitions``.

       In general, the more administrative work you can move into your functions
       the better.  That way the Dask scheduler doesn't need to think about all
       of the fine-grained operations.

    -  **Breaking up your computation:** For very large workloads you may also want to
       try sending smaller chunks to Dask at a time.  For example if you're
       processing a petabyte of data but find that Dask is only happy with 100
       TB, maybe you can break up your computation into ten pieces and submit
       them one after the other.


Learn Techniques For Customization
----------------------------------

The high level Dask collections (array, dataframe, bag) include common
operations that follow standard Python APIs from NumPy and Pandas.
However, many Python workloads are complex and may require operations that are
not included in these high level APIs.

Fortunately, there are many options to support custom workloads:

-   All collections have a ``map_partitions`` or ``map_blocks`` function, that
    applies a user provided function across every Pandas dataframe or NumPy array
    in the collection.  Because Dask collections are made up of normal Python
    objects, it's often quite easy to map custom functions across partitions of a
    dataset without much modification.

    .. code-block:: python

       df.map_partitions(my_custom_func)

-   More complex ``map_*`` functions.  Sometimes your custom behavior isn't
    embarrassingly parallel, but requires more advanced communication.  For
    example maybe you need to communicate a little bit of information from one
    partition to the next, or maybe you want to build a custom aggregation.

    Dask collections include methods for these as well.

-   For even more complex workloads you can convert your collections into
    individual blocks, and arrange those blocks as you like using Dask Delayed.
    There is usually a ``to_delayed`` method on every collection.

.. currentmodule:: dask.dataframe

.. autosummary::

    map_partitions
    rolling.map_overlap
    groupby.Aggregation

.. currentmodule:: dask.array

.. autosummary::

    blockwise
    map_blocks
    map_overlap
    reduction


Stop Using Dask When No Longer Needed
-------------------------------------

In many workloads it is common to use Dask to read in a large amount of data,
reduce it down, and then iterate on a much smaller amount of data.  For this
latter stage on smaller data it may make sense to stop using Dask, and start
using normal Python again.

.. code-block:: python

   df = dd.read_parquet("lots-of-data-*.parquet")
   df = df.groupby('name').mean()  # reduce data significantly
   df = df.compute()               # continue on with Pandas/NumPy


Persist When You Can
--------------------

Accessing data from RAM is often much faster than accessing it from disk.
Once you have your dataset in a clean state that both:

1.  Fits in memory
2.  Is clean enough that you will want to try many different analyses

Then it is a good time to *persist* your data in RAM

.. code-block:: python

    df = dd.read_parquet("lots-of-data-*.parquet")
    df = df.fillna(...)  # clean up things lazily
    df = df[df.name == 'Alice']  # get down to a more reasonable size

    df = df.persist()  # trigger computation, persist in distributed RAM

Note that this is only relevant if you are on a distributed machine (otherwise,
as mentioned above, you should probably continue on without Dask).


Store Data Efficiently
----------------------

As your ability to compute increases you will likely find that data access and
I/O take up a larger portion of your total time.  Additionally, parallel
computing will often add new constraints to how your store your data,
particularly around providing random access to blocks of your data that are in
line with how you plan to compute on it.

For example ...

-   For compression you'll probably find that you drop gzip and bz2, and embrace
    newer systems like lz4, snappy, and Z-Standard that provide better
    performance and random access.
-   For storage formats you may find that you want self-describing formats that
    are optimized for random access, metadata storage, and binary encoding like
    Parquet, ORC, Zarr, HDF5, GeoTIFF and so on
-   When working on the cloud you may find that some older formats like HDF5 may
    not work well
-   You may want to partition or chunk your data in ways that align well to
    common queries.  In Dask DataFrame this might mean choosing a column to
    sort by for fast selection and joins.  For Dask Array this might mean
    choosing chunk sizes that are aligned with your access patterns and
    algorithms.

Processes and Threads
---------------------

If you're doing mostly numeric work with Numpy, Pandas, Scikit-Learn, Numba,
and other libraries that release the `GIL <https://docs.python.org/3/glossary.html#term-global-interpreter-lock>`_, then use mostly threads.  If you're
doing work on text data or Python collections like lists and dicts then use
mostly processes.

If you're on larger machines with a high thread count (greater than 10), then
you should probably split things up into at least a few processes regardless.
Python can be highly productive with 10 threads per process with numeric work,
but not 50 threads.

For more information on threads, processes, and how to configure them in Dask, see
:doc:`the scheduler documentation <scheduling>`.


Load Data with Dask
-------------------

If you need to work with large Python objects, then please let Dask create
them.  A common anti-pattern we see is people creating large Python objects
outside of Dask, then giving those objects to Dask and asking it to manage them.
This works, but means that Dask needs to move around these very large objects
with its metadata, rather than as normal Dask-controlled results.

Here are some common patterns to avoid and nicer alternatives:

DataFrames
~~~~~~~~~~

.. code-block:: python

   # Don't

   ddf = ... a dask dataframe ...
   for fn in filenames:
       df = pandas.read_csv(fn)  # Read locally with Pandas
       ddf = ddf.append(df)            # Give to Dask

.. code-block:: python

    # Do

    ddf = dd.read_csv(filenames)

Arrays
~~~~~~

.. code-block:: python

   # Don't

   f = h5py.File(...)
   x = np.asarray(f["x"])  # Get data as a NumPy array locally

   x = da.from_array(x)  # Hand NumPy array to Dask

.. code-block:: python

   # Do

   f = h5py.File(...)
   x = da.from_array(f["x"])  # Let Dask do the reading

Delayed
~~~~~~~

.. code-block:: python

    # Don't

    @dask.delayed
    def process(a, b):
        ...

    df = pandas.read_csv("some-large-file.csv")  # Create large object locally
    results = []
    for item in L:
        result = process(item, df)  # include df in every delayed call
        results.append(result)

.. code-block:: python

   # Do

   @dask.delayed
   def process(a, b):
       ...

   df = dask.delayed(pandas.read_csv)("some-large-file.csv")  # Let Dask build object
   results = []
   for item in L:
       result = process(item, df)  # include pointer to df in every delayed call
       results.append(result)


Avoid calling compute repeatedly
--------------------------------

Compute related results with shared computations in a single :func:`dask.compute` call

.. code-block:: python

   # Don't repeatedly call compute

   df = dd.read_csv("...")
   xmin = df.x.min().compute()
   xmax = df.x.max().compute()

.. code-block:: python

   # Do compute multiple results at the same time

   df = dd.read_csv("...")

   xmin, xmax = dask.compute(df.x.min(), df.x.max())

This allows Dask to compute the shared parts of the computation (like the
``dd.read_csv`` call above) only once, rather than once per ``compute`` call.
High Level Graphs
=================

Dask graphs produced by collections like Arrays, Bags, and DataFrames have
high-level structure that can be useful for visualization and high-level
optimization.  The task graphs produced by these collections encode this
structure explicitly as ``HighLevelGraph`` objects.  This document describes
how to work with these in more detail.


Motivation and Example
----------------------

In full generality, Dask schedulers expect arbitrary task graphs where each
node is a single Python function call and each edge is a dependency between
two function calls.  These are usually stored in flat dictionaries.  Here is
some simple Dask DataFrame code and the task graph that it might generate:

.. code-block:: python

    import dask.dataframe as dd

    df = dd.read_csv('myfile.*.csv')
    df = df + 100
    df = df[df.name == 'Alice']

.. code-block:: python

   {
    ('read-csv', 0): (pandas.read_csv, 'myfile.0.csv'),
    ('read-csv', 1): (pandas.read_csv, 'myfile.1.csv'),
    ('read-csv', 2): (pandas.read_csv, 'myfile.2.csv'),
    ('read-csv', 3): (pandas.read_csv, 'myfile.3.csv'),
    ('add', 0): (operator.add, ('read-csv', 0), 100),
    ('add', 1): (operator.add, ('read-csv', 1), 100),
    ('add', 2): (operator.add, ('read-csv', 2), 100),
    ('add', 3): (operator.add, ('read-csv', 3), 100),
    ('filter', 0): (lambda part: part[part.name == 'Alice'], ('add', 0)),
    ('filter', 1): (lambda part: part[part.name == 'Alice'], ('add', 1)),
    ('filter', 2): (lambda part: part[part.name == 'Alice'], ('add', 2)),
    ('filter', 3): (lambda part: part[part.name == 'Alice'], ('add', 3)),
   }

The task graph is a dictionary that stores every Pandas-level function call
necessary to compute the final result.  We can see that there is some structure
to this dictionary if we separate out the tasks that were associated to each
high-level Dask DataFrame operation:

.. code-block:: python

   {
    # From the dask.dataframe.read_csv call
    ('read-csv', 0): (pandas.read_csv, 'myfile.0.csv'),
    ('read-csv', 1): (pandas.read_csv, 'myfile.1.csv'),
    ('read-csv', 2): (pandas.read_csv, 'myfile.2.csv'),
    ('read-csv', 3): (pandas.read_csv, 'myfile.3.csv'),

    # From the df + 100 call
    ('add', 0): (operator.add, ('read-csv', 0), 100),
    ('add', 1): (operator.add, ('read-csv', 1), 100),
    ('add', 2): (operator.add, ('read-csv', 2), 100),
    ('add', 3): (operator.add, ('read-csv', 3), 100),

    # From the df[df.name == 'Alice'] call
    ('filter', 0): (lambda part: part[part.name == 'Alice'], ('add', 0)),
    ('filter', 1): (lambda part: part[part.name == 'Alice'], ('add', 1)),
    ('filter', 2): (lambda part: part[part.name == 'Alice'], ('add', 2)),
    ('filter', 3): (lambda part: part[part.name == 'Alice'], ('add', 3)),
   }

By understanding this high-level structure we are able to understand our task
graphs more easily (this is more important for larger datasets when there are
thousands of tasks per layer) and how to perform high-level optimizations.  For
example, in the case above we may want to automatically rewrite our code to
filter our datasets before adding 100:

.. code-block:: python

    # Before
    df = dd.read_csv('myfile.*.csv')
    df = df + 100
    df = df[df.name == 'Alice']

    # After
    df = dd.read_csv('myfile.*.csv')
    df = df[df.name == 'Alice']
    df = df + 100

Dask's high level graphs help us to explicitly encode this structure by storing
our task graphs in layers with dependencies between layers:

.. code-block:: python

   >>> import dask.dataframe as dd

   >>> df = dd.read_csv('myfile.*.csv')
   >>> df = df + 100
   >>> df = df[df.name == 'Alice']

   >>> graph = df.__dask_graph__()
   >>> graph.layers
   {
    'read-csv': {('read-csv', 0): (pandas.read_csv, 'myfile.0.csv'),
                 ('read-csv', 1): (pandas.read_csv, 'myfile.1.csv'),
                 ('read-csv', 2): (pandas.read_csv, 'myfile.2.csv'),
                 ('read-csv', 3): (pandas.read_csv, 'myfile.3.csv')},

    'add': {('add', 0): (operator.add, ('read-csv', 0), 100),
            ('add', 1): (operator.add, ('read-csv', 1), 100),
            ('add', 2): (operator.add, ('read-csv', 2), 100),
            ('add', 3): (operator.add, ('read-csv', 3), 100)}

    'filter': {('filter', 0): (lambda part: part[part.name == 'Alice'], ('add', 0)),
               ('filter', 1): (lambda part: part[part.name == 'Alice'], ('add', 1)),
               ('filter', 2): (lambda part: part[part.name == 'Alice'], ('add', 2)),
               ('filter', 3): (lambda part: part[part.name == 'Alice'], ('add', 3))}
   }

   >>> graph.dependencies
   {
    'read-csv': set(),
    'add': {'read-csv'},
    'filter': {'add'}
   }

While the DataFrame points to the output layers on which it depends directly:

.. code-block:: python

   >>> df.__dask_layers__()
   {'filter'}


HighLevelGraphs
---------------

The :obj:`HighLevelGraph` object is a ``Mapping`` object composed of other
sub-``Mappings``, along with a high-level dependency mapping between them:

.. code-block:: python

   class HighLevelGraph(Mapping):
       layers: Dict[str, Mapping]
       dependencies: Dict[str, Set[str]]

You can construct a HighLevelGraph explicitly by providing both to the
constructor:

.. code-block:: python

   layers = {
      'read-csv': {('read-csv', 0): (pandas.read_csv, 'myfile.0.csv'),
                   ('read-csv', 1): (pandas.read_csv, 'myfile.1.csv'),
                   ('read-csv', 2): (pandas.read_csv, 'myfile.2.csv'),
                   ('read-csv', 3): (pandas.read_csv, 'myfile.3.csv')},

      'add': {('add', 0): (operator.add, ('read-csv', 0), 100),
              ('add', 1): (operator.add, ('read-csv', 1), 100),
              ('add', 2): (operator.add, ('read-csv', 2), 100),
              ('add', 3): (operator.add, ('read-csv', 3), 100)},

      'filter': {('filter', 0): (lambda part: part[part.name == 'Alice'], ('add', 0)),
                 ('filter', 1): (lambda part: part[part.name == 'Alice'], ('add', 1)),
                 ('filter', 2): (lambda part: part[part.name == 'Alice'], ('add', 2)),
                 ('filter', 3): (lambda part: part[part.name == 'Alice'], ('add', 3))}
   }

   dependencies = {'read-csv': set(),
                   'add': {'read-csv'},
                   'filter': {'add'}}

   graph = HighLevelGraph(layers, dependencies)

This object satisfies the ``Mapping`` interface, and so operates as a normal
Python dictionary that is the semantic merger of the underlying layers:

.. code-block:: python

   >>> len(graph)
   12
   >>> graph[('read-csv', 0)]
   ('read-csv', 0): (pandas.read_csv, 'myfile.0.csv'),


API
---

.. currentmodule:: dask.highlevelgraph

.. autoclass:: HighLevelGraph
   :members:
   :inherited-members:
Create Dask Arrays
==================

You can load or store Dask arrays from a variety of common sources like HDF5,
NetCDF, `Zarr`_, or any format that supports NumPy-style slicing.

.. currentmodule:: dask.array

.. autosummary::
   from_array
   from_delayed
   from_npy_stack
   from_zarr
   stack
   concatenate

NumPy Slicing
-------------

.. autosummary::
   from_array

Many storage formats have Python projects that expose storage using NumPy
slicing syntax.  These include HDF5, NetCDF, BColz, Zarr, GRIB, etc.  For
example, we can load a Dask array from an HDF5 file using `h5py <https://www.h5py.org/>`_:

.. code-block:: Python

   >>> import h5py
   >>> f = h5py.File('myfile.hdf5') # HDF5 file
   >>> d = f['/data/path']          # Pointer on on-disk array
   >>> d.shape                      # d can be very large
   (1000000, 1000000)

   >>> x = d[:5, :5]                # We slice to get numpy arrays

Given an object like ``d`` above that has ``dtype`` and ``shape`` properties
and that supports NumPy style slicing, we can construct a lazy Dask array:

.. code-block:: Python

   >>> import dask.array as da
   >>> x = da.from_array(d, chunks=(1000, 1000))

This process is entirely lazy.  Neither creating the h5py object nor wrapping
it with ``da.from_array`` have loaded any data.


Random Data
-----------

For experimentation or benchmarking it is common to create arrays of random
data.  The ``dask.array.random`` module implements most of the functions in the
``numpy.random`` module.  We list some common functions below but for a full
list see the :doc:`Array API <array-api>`:

.. autosummary::
   random.binomial
   random.normal
   random.poisson
   random.random

.. code-block:: python

   >>> import dask.array as da
   >>> x = da.random.random((10000, 10000), chunks=(1000, 1000))


Concatenation and Stacking
--------------------------

.. autosummary::
   stack
   concatenate

Often we store data in several different locations and want to stitch them together:

.. code-block:: Python

    dask_arrays = []
    for fn in filenames:
        f = h5py.File(fn)
        d = f['/data']
        array = da.from_array(d, chunks=(1000, 1000))
        dask_arrays.append(array)

    x = da.concatenate(dask_arrays, axis=0)  # concatenate arrays along first axis

For more information, see :doc:`concatenation and stacking <array-stack>` docs.


Using ``dask.delayed``
----------------------

.. autosummary::
   from_delayed
   stack
   concatenate

Sometimes NumPy-style data resides in formats that do not support NumPy-style
slicing.  We can still construct Dask arrays around this data if we have a
Python function that can generate pieces of the full array if we use
:doc:`dask.delayed <delayed>`.  Dask delayed lets us delay a single function
call that would create a NumPy array.  We can then wrap this delayed object
with ``da.from_delayed``, providing a dtype and shape to produce a
single-chunked Dask array.  Furthermore, we can use ``stack`` or ``concatenate`` from
before to construct a larger lazy array.

As an example, consider loading a stack of images using ``skimage.io.imread``:

.. code-block:: python

    import skimage.io
    import dask.array as da
    import dask

    imread = dask.delayed(skimage.io.imread, pure=True)  # Lazy version of imread

    filenames = sorted(glob.glob('*.jpg'))

    lazy_images = [imread(path) for path in filenames]   # Lazily evaluate imread on each path
    sample = lazy_images[0].compute()  # load the first image (assume rest are same shape/dtype)

    arrays = [da.from_delayed(lazy_image,           # Construct a small Dask array
                              dtype=sample.dtype,   # for every lazy value
                              shape=sample.shape)
              for lazy_image in lazy_images]

    stack = da.stack(arrays, axis=0)                # Stack all small Dask arrays into one

See :doc:`documentation on using dask.delayed with collections<delayed-collections>`.

Often it is substantially faster to use ``da.map_blocks`` rather than ``da.stack``

.. code-block:: python

    import glob
    import skimage.io
    import numpy as np
    import dask.array as da

    filenames = sorted(glob.glob('*.jpg'))

    def read_one_image(block_id, filenames=filenames, axis=0):
        # a function that reads in one chunk of data
        path = filenames[block_id[axis]]
        image = skimage.io.imread(path)
        return np.expand_dims(image, axis=axis)

    # load the first image (assume rest are same shape/dtype)
    sample = skimage.io.imread(filenames[0])

    stack = da.map_blocks(
        read_one_image,
        dtype=sample.dtype,
        chunks=((1,) * len(filenames),  *sample.shape)
    )


From Dask DataFrame
-------------------

There are several ways to create a Dask array from a Dask DataFrame. Dask DataFrames have a ``to_dask_array`` method:

.. code-block:: python

   >>> df = dask.dataframes.from_pandas(...)
   >>> df.to_dask_array()
   dask.array<values, shape=(nan, 3), dtype=float64, chunksize=(nan, 3), chunktype=numpy.ndarray>

This mirrors the `to_numpy
<https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_numpy.html>`_
function in Pandas. The ``values`` attribute is also supported:

.. code-block:: python

   >>> df.values
   dask.array<values, shape=(nan, 3), dtype=float64, chunksize=(nan, 3), chunktype=numpy.ndarray>

However, these arrays do not have known chunk sizes because dask.dataframe does
not track the number of rows in each partition. This means that some operations
like slicing will not operate correctly.

The chunk sizes can be computed:

.. code-block:: python

   >>> df.to_dask_array(lengths=True)
   dask.array<array, shape=(100, 3), dtype=float64, chunksize=(50, 3), chunktype=numpy.ndarray>

Specifying ``lengths=True`` triggers immediate computation of the chunk sizes.
This enables downstream computations that rely on having known chunk sizes
(e.g., slicing).

The Dask DataFrame ``to_records`` method also returns a Dask Array, but does not compute the shape
information:

.. code-block:: python

   >>> df.to_records()
   dask.array<to_records, shape=(nan,), dtype=(numpy.record, [('index', '<i8'), ('x', '<f8'), ('y', '<f8'), ('z', '<f8')]), chunksize=(nan,), chunktype=numpy.ndarray>

If you have a function that converts a Pandas DataFrame into a NumPy array,
then calling ``map_partitions`` with that function on a Dask DataFrame will
produce a Dask array:

.. code-block:: python

   >>> df.map_partitions(np.asarray)
   dask.array<asarray, shape=(nan, 3), dtype=float64, chunksize=(nan, 3), chunktype=numpy.ndarray>


Interactions with NumPy arrays
------------------------------

Dask array operations will automatically convert NumPy arrays into single-chunk
dask arrays:

.. code-block:: python

   >>> x = da.sum(np.ones(5))
   >>> x.compute()
   5

When NumPy and Dask arrays interact, the result will be a Dask array.  Automatic
rechunking rules will generally slice the NumPy array into the appropriate Dask
chunk shape:

.. code-block:: python

   >>> x = da.ones(10, chunks=(5,))
   >>> y = np.ones(10)
   >>> z = x + y
   >>> z
   dask.array<add, shape=(10,), dtype=float64, chunksize=(5,), chunktype=numpy.ndarray>

These interactions work not just for NumPy arrays but for any object that has
shape and dtype attributes and implements NumPy slicing syntax.

Memory mapping
--------------

Memory mapping can be a highly effective method to access raw binary data since
it has nearly zero overhead if the data is already in the file system cache. For
the threaded scheduler, creating a Dask array from a raw binary file can be as simple as
:code:`a = da.from_array(np.memmap(filename, shape=shape, dtype=dtype, mode='r'))`.

For multiprocessing or distributed schedulers, the memory map for each array
chunk should be created on the correct worker process and not on the main
process to avoid data transfer through the cluster. This can be achieved by
wrapping the function that creates the memory map using :code:`dask.delayed`.

.. code-block:: python

   import numpy as np
   import dask
   import dask.array as da


   def mmap_load_chunk(filename, shape, dtype, offset, sl):
       '''
       Memory map the given file with overall shape and dtype and return a slice
       specified by :code:`sl`.

       Parameters
       ----------

       filename : str
       shape : tuple
           Total shape of the data in the file
       dtype:
           NumPy dtype of the data in the file
       offset : int
           Skip :code:`offset` bytes from the beginning of the file.
       sl:
           Object that can be used for indexing or slicing a NumPy array to
           extract a chunk

       Returns
       -------

       numpy.memmap or numpy.ndarray
           View into memory map created by indexing with :code:`sl`,
           or NumPy ndarray in case no view can be created using :code:`sl`.
       '''
       data = np.memmap(filename, mode='r', shape=shape, dtype=dtype, offset=offset)
       return data[sl]


   def mmap_dask_array(filename, shape, dtype, offset=0, blocksize=5):
       '''
       Create a Dask array from raw binary data in :code:`filename`
       by memory mapping.

       This method is particularly effective if the file is already
       in the file system cache and if arbitrary smaller subsets are
       to be extracted from the Dask array without optimizing its
       chunking scheme.

       It may perform poorly on Windows if the file is not in the file
       system cache. On Linux it performs well under most circumstances.

       Parameters
       ----------

       filename : str
       shape : tuple
           Total shape of the data in the file
       dtype:
           NumPy dtype of the data in the file
       offset : int, optional
           Skip :code:`offset` bytes from the beginning of the file.
       blocksize : int, optional
           Chunk size for the outermost axis. The other axes remain unchunked.

       Returns
       -------

       dask.array.Array
           Dask array matching :code:`shape` and :code:`dtype`, backed by
           memory-mapped chunks.
       '''
       load = dask.delayed(mmap_load_chunk)
       chunks = []
       for index in range(0, shape[0], blocksize):
           # Truncate the last chunk if necessary
           chunk_size = min(blocksize, shape[0] - index)
           chunk = dask.array.from_delayed(
               load(
                   filename,
                   shape=shape,
                   dtype=dtype,
                   offset=offset,
                   sl=slice(index, index + chunk_size)
               ),
               shape=(chunk_size, ) + shape[1:],
               dtype=dtype
           )
           chunks.append(chunk)
       return da.concatenate(chunks, axis=0)

   x = mmap_dask_array(
       filename='testfile-50-50-100-100-float32.raw',
       shape=(50, 50, 100, 100),
       dtype=np.float32
   )


Chunks
------

See :doc:`documentation on Array Chunks <array-chunks>` for more information.


Store Dask Arrays
=================

.. autosummary::
   store
   to_hdf5
   to_npy_stack
   to_zarr
   compute

In Memory
---------

.. autosummary::
   compute

If you have a small amount of data, you can call ``np.array`` or ``.compute()``
on your Dask array to turn in to a normal NumPy array:

.. code-block:: Python

   >>> x = da.arange(6, chunks=3)
   >>> y = x**2
   >>> np.array(y)
   array([0, 1, 4, 9, 16, 25])

   >>> y.compute()
   array([0, 1, 4, 9, 16, 25])


NumPy style slicing
-------------------

.. autosummary::
   store

You can store Dask arrays in any object that supports NumPy-style slice
assignment like ``h5py.Dataset``:

.. code-block:: Python

   >>> import h5py
   >>> f = h5py.File('myfile.hdf5')
   >>> d = f.require_dataset('/data', shape=x.shape, dtype=x.dtype)
   >>> da.store(x, d)

Also, you can store several arrays in one computation by passing lists of sources and
destinations:

.. code-block:: Python

   >>> da.store([array1, array2], [output1, output2])  # doctest: +SKIP

HDF5
----

.. autosummary::
   to_hdf5

HDF5 is sufficiently common that there is a special function ``to_hdf5`` to
store data into HDF5 files using ``h5py``:

.. code-block:: Python

   >>> da.to_hdf5('myfile.hdf5', '/y', y)  # doctest: +SKIP

You can store several arrays in one computation with the function
``da.to_hdf5`` by passing in a dictionary:

.. code-block:: Python

   >>> da.to_hdf5('myfile.hdf5', {'/x': x, '/y': y})  # doctest: +SKIP


Zarr
----

The `Zarr`_ format is a chunk-wise binary array
storage file format with a good selection of encoding and compression options.
Due to each chunk being stored in a separate file, it is ideal for parallel
access in both reading and writing (for the latter, if the Dask array
chunks are aligned with the target). Furthermore, storage in
:doc:`remote data services <how-to/connect-to-remote-data>` such as S3 and GCS is
supported.

For example, to save data to a local zarr dataset you would do:

.. code-block:: Python

   >>> arr.to_zarr('output.zarr')

or to save to a particular bucket on S3:

.. code-block:: Python

   >>> arr.to_zarr('s3://mybucket/output.zarr', storage_option={'key': 'mykey',
                   'secret': 'mysecret'})

or your own custom zarr Array:

.. code-block:: Python

   >>> z = zarr.create((10,), dtype=float, store=zarr.ZipStore("output.zarr"))
   >>> arr.to_zarr(z)

To retrieve those data, you would do ``da.from_zarr`` with exactly the same arguments. The
chunking of the resultant Dask array is defined by how the files were saved, unless
otherwise specified.


TileDB
------

`TileDB <https://docs.tiledb.io>`_  is a binary array format and storage manager with
tunable chunking, layout, and compression options. The TileDB storage manager library
includes support for scalable storage backends such as S3 API compatible object stores
and HDFS, with automatic scaling, and supports multi-threaded and multi-process
reads (consistent) and writes (eventually-consistent).

To save data to a local TileDB array:

.. code-block:: Python

  >>> arr.to_tiledb('output.tdb')

or to save to a bucket on S3:

.. code-block:: python

  >>> arr.to_tiledb('s3://mybucket/output.tdb',
                    storage_options={'vfs.s3.aws_access_key_id': 'mykey',
                                     'vfs.s3.aws_secret_access_key': 'mysecret'})

Files may be retrieved by running `da.from_tiledb` with the same URI, and any
necessary arguments.


Intermediate storage
--------------------

.. autosummary::
   store

In some cases, one may wish to store an intermediate result in long term
storage. This differs from ``persist``, which is mainly used to manage
intermediate results within Dask that don't necessarily have longevity.
Also it differs from storing final results as these mark the end of the Dask
graph. Thus intermediate results are easier to reuse without reloading data.
Intermediate storage is mainly useful in cases where the data is needed
outside of Dask (e.g. on disk, in a database, in the cloud, etc.). It can
be useful as a checkpoint for long running or error-prone computations.

The intermediate storage use case differs from the typical storage use case as
a Dask Array is returned to the user that represents the result of that
storage operation. This is typically done by setting the ``store`` function's
``return_stored`` flag to ``True``.

.. code-block:: python

   x.store()  # stores data, returns nothing
   x = x.store(return_stored=True)  # stores data, returns new dask array backed by that data

The user can then decide whether the
storage operation happens immediately (by setting the ``compute`` flag to
``True``) or later (by setting the ``compute`` flag to ``False``). In all
other ways, this behaves the same as a normal call to ``store``. Some examples
are shown below.

.. code-block:: Python

   >>> import dask.array as da
   >>> import zarr as zr
   >>> c = (2, 2)
   >>> d = da.ones((10, 11), chunks=c)
   >>> z1 = zr.open_array('lazy.zarr', shape=d.shape, dtype=d.dtype, chunks=c)
   >>> z2 = zr.open_array('eager.zarr', shape=d.shape, dtype=d.dtype, chunks=c)
   >>> d1 = d.store(z1, compute=False, return_stored=True)
   >>> d2 = d.store(z2, compute=True, return_stored=True)

This can be combined with any other storage strategies either noted above, in
the docs or for any specialized storage types.


Plugins
=======

We can run arbitrary user-defined functions on Dask arrays whenever they are
constructed.  This allows us to build a variety of custom behaviors that improve
debugging, user warning, etc.  You can register a list of functions to run on
all Dask arrays to the global ``array_plugins=`` value:

.. code-block:: python

   >>> def f(x):
   ...     print(x.nbytes)

   >>> with dask.config.set(array_plugins=[f]):
   ...     x = da.ones((10, 1), chunks=(5, 1))
   ...     y = x.dot(x.T)
   80
   80
   800
   800

If the plugin function returns None, then the input Dask array will be returned
without change.  If the plugin function returns something else, then that value
will be the result of the constructor.

Examples
--------

Automatically compute
~~~~~~~~~~~~~~~~~~~~~

We may wish to turn some Dask array code into normal NumPy code.  This is
useful, for example, to track down errors immediately that would otherwise be
hidden by Dask's lazy semantics:

.. code-block:: python

   >>> with dask.config.set(array_plugins=[lambda x: x.compute()]):
   ...     x = da.arange(5, chunks=2)

   >>> x  # this was automatically converted into a numpy array
   array([0, 1, 2, 3, 4])

Warn on large chunks
~~~~~~~~~~~~~~~~~~~~

We may wish to warn users if they are creating chunks that are too large:

.. code-block:: python

   def warn_on_large_chunks(x):
       shapes = list(itertools.product(*x.chunks))
       nbytes = [x.dtype.itemsize * np.prod(shape) for shape in shapes]
       if any(nb > 1e9 for nb in nbytes):
           warnings.warn("Array contains very large chunks")

   with dask.config.set(array_plugins=[warn_on_large_chunks]):
       ...

Combine
~~~~~~~

You can also combine these plugins into a list.  They will run one after the
other, chaining results through them:

.. code-block:: python

   with dask.config.set(array_plugins=[warn_on_large_chunks, lambda x: x.compute()]):
       ...


.. _Zarr: https://zarr.readthedocs.io/en/stable/
:orphan:

Scheduler Overview
==================

After we create a dask graph, we use a scheduler to run it. Dask currently
implements a few different schedulers:

-  ``dask.threaded.get``: a scheduler backed by a thread pool
-  ``dask.multiprocessing.get``: a scheduler backed by a process pool
-  ``dask.get``: a synchronous scheduler, good for debugging
-  ``distributed.Client.get``: a distributed scheduler for executing graphs
   on multiple machines.  This lives in the external distributed_ project.

.. _distributed: https://distributed.dask.org/en/latest/


The ``get`` function
--------------------

The entry point for all schedulers is a ``get`` function. This takes a dask
graph, and a key or list of keys to compute:

.. code-block:: python

   >>> from operator import add

   >>> dsk = {'a': 1,
   ...        'b': 2,
   ...        'c': (add, 'a', 'b'),
   ...        'd': (sum, ['a', 'b', 'c'])}

   >>> get(dsk, 'c')
   3

   >>> get(dsk, 'd')
   6

   >>> get(dsk, ['a', 'b', 'c'])
   [1, 2, 3]


Using ``compute`` methods
-------------------------

When working with dask collections, you will rarely need to
interact with scheduler ``get`` functions directly. Each collection has a
default scheduler, and a built-in ``compute`` method that calculates the output
of the collection:

.. code-block:: python

    >>> import dask.array as da
    >>> x = da.arange(100, chunks=10)
    >>> x.sum().compute()
    4950

The compute method takes a number of keywords:

- ``scheduler``: the name of the desired scheduler as a string (``"threads"``, ``"processes"``, ``"single-threaded"``, etc.), a ``get`` function, or a ``dask.distributed.Client`` object.  Overrides the default for the collection.
- ``**kwargs``: extra keywords to pass on to the scheduler ``get`` function.

See also: :ref:`configuring-schedulers`.


The ``compute`` function
------------------------

You may wish to compute results from multiple dask collections at once.
Similar to the ``compute`` method on each collection, there is a general
``compute`` function that takes multiple collections and returns multiple
results. This merges the graphs from each collection, so intermediate results
are shared:

.. code-block:: python

    >>> y = (x + 1).sum()
    >>> z = (x + 1).mean()
    >>> da.compute(y, z)    # Compute y and z, sharing intermediate results
    (5050, 50.5)

Here the ``x + 1`` intermediate was only computed once, while calling
``y.compute()`` and ``z.compute()`` would compute it twice. For large graphs
that share many intermediates, this can be a big performance gain.

The ``compute`` function works with any dask collection, and is found in
``dask.base``. For convenience it has also been imported into the top level
namespace of each collection.

.. code-block:: python

    >>> from dask.base import compute
    >>> compute is da.compute
    True


.. _configuring-schedulers:

Configuring the schedulers
--------------------------

The dask collections each have a default scheduler:

- ``dask.array`` and ``dask.dataframe`` use the threaded scheduler by default
- ``dask.bag`` uses the multiprocessing scheduler by default.

For most cases, the default settings are good choices. However, sometimes you
may want to use a different scheduler. There are two ways to do this.

1. Using the ``scheduler`` keyword in the ``compute`` method:

    .. code-block:: python

        >>> x.sum().compute(scheduler='processes')

2. Using ``dask.config.set``. This can be used either as a context manager, or to
   set the scheduler globally:

    .. code-block:: python

        # As a context manager
        >>> with dask.config.set(scheduler='processes'):
        ...     x.sum().compute()

        # Set globally
        >>> dask.config.set(scheduler='processes')
        >>> x.sum().compute()


Additionally, each scheduler may take a few extra keywords specific to that
scheduler. For example, the multiprocessing and threaded schedulers each take a
``num_workers`` keyword, which sets the number of processes or threads to use
(defaults to number of cores). This can be set by passing the keyword when
calling ``compute``:

.. code-block:: python

    # Compute with 4 threads
    >>> x.compute(num_workers=4)

Alternatively, the multiprocessing and threaded schedulers will check for a
global pool set with ``dask.config.set``:

.. code-block:: python

    >>> from concurrent.futures import ThreadPoolExecutor
    >>> with dask.config.set(pool=ThreadPoolExecutor(4)):
    ...     x.compute()

The multiprocessing scheduler also supports `different contexts`_ ("spawn",
"forkserver", "fork") which you can set with ``dask.config.set``. The default
context is "spawn", but you can set a different one:

.. code-block:: python

   >>> with dask.config.set({"multiprocessing.context": "forkserver"}):
   ...     x.compute()

.. _different contexts: https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods

For more information on the individual options for each scheduler, see the
docstrings for each scheduler ``get`` function.


Debugging the schedulers
------------------------

Debugging parallel code can be difficult, as conventional tools such as ``pdb``
don't work well with multiple threads or processes. To get around this when
debugging, we recommend using the synchronous scheduler found at
``dask.get``. This runs everything serially, allowing it to work
well with ``pdb``:

.. code-block:: python

    >>> dask.config.set(scheduler='single-threaded')
    >>> x.sum().compute()    # This computation runs serially instead of in parallel


The shared memory schedulers also provide a set of callbacks that can be used
for diagnosing and profiling. You can learn more about scheduler callbacks and
diagnostics :doc:`here <diagnostics-local>`.


More Information
----------------

- See :doc:`shared` for information on the design of the shared memory
  (threaded or multiprocessing) schedulers
- See distributed_ for information on the distributed memory scheduler
FAQ
===

**Question**: *Is Dask appropriate for adoption within a larger institutional context?*

**Answer**: *Yes.* Dask is used within the world's largest banks, national labs,
retailers, technology companies, and government agencies.  It is used in highly
secure environments.  It is used in conservative institutions as well as fast
moving ones.

This page contains Frequently Asked Questions and concerns from institutions and users
when they first investigate Dask.

.. contents:: :local:

For Management
--------------

Briefly, what problem does Dask solve for us?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask is a general purpose parallel programming solution.
As such it is used in *many* different ways.

However, the most common problem that Dask solves is connecting Python analysts
to distributed hardware, particularly for data science and machine learning
workloads.  The institutions for whom Dask has the greatest
impact are those who have a large body of Python users who are accustomed to
libraries like NumPy, Pandas, Jupyter, Scikit-Learn and others, but want to
scale those workloads across a cluster.  Often they also have distributed
computing resources that are going underused.

Dask removes both technological and cultural barriers to connect Python users
to computing resources in a way that is native to both the users and IT.

"*Help me scale my notebook onto the cluster*" is a common pain point for
institutions today, and it is a common entry point for Dask usage.


Is Dask mature?  Why should we trust it?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes.  While Dask itself is relatively new (it began in 2015) it is built by the
NumPy, Pandas, Jupyter, Scikit-Learn developer community, which is well trusted.
Dask is a relatively thin wrapper on top of these libraries and,
as a result, the project can be relatively small and simple.
It doesn't reinvent a whole new system.

Additionally, this tight integration with the broader technology stack
gives substantial benefits long term.  For example:

-   Because Pandas maintainers also maintain Dask,
    when Pandas issues a new releases Dask issues a release at the same time to
    ensure continuity and compatibility.

-   Because Scikit-Learn maintainers maintain and use Dask when they train on large clusters,
    you can be assured that Dask-ML focuses on pragmatic and important
    solutions like XGBoost integration, and hyper-parameter selection,
    and that the integration between the two feels natural for novice and
    expert users alike.

-   Because Jupyter maintainers also maintain Dask,
    powerful Jupyter technologies like JupyterHub and JupyterLab are designed
    with Dask's needs in mind, and new features are pushed quickly to provide a
    first class and modern user experience.

Additionally, Dask is maintained both by a broad community of maintainers,
as well as substantial institutional support (several full-time employees each)
by both Anaconda, the company behind the leading data science distribution, and
NVIDIA, the leading hardware manufacturer of GPUs.  Despite large corporate
support, Dask remains a community governed project, and is fiscally sponsored by
NumFOCUS, the same 501c3 that fiscally sponsors NumPy, Pandas, Jupyter, and many others.


Who else uses Dask?
~~~~~~~~~~~~~~~~~~~

Dask is used by individual researchers in practically every field today.  It
has millions of downloads per month, and is integrated into many PyData
software packages today.

On an *institutional* level Dask is used by analytics and research groups in a
similarly broad set of domains across both energetic startups as well as large
conservative household names.  A web search shows articles by Capital One,
Barclays, Walmart, NASA, Los Alamos National Laboratories, and hundreds of
other similar institutions.


How does Dask compare with Apache Spark?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*This question has longer and more technical coverage* :doc:`here <spark>`

Dask and Apache Spark are similar in that they both ...

-  Promise easy parallelism for data science Python users
-  Provide Dataframe and ML APIs for ETL, data science, and machine learning
-  Scale out to similar scales, around 1-1000 machines

Dask differs from Apache Spark in a few ways:

-  Dask is more Python native, Spark is Scala/JVM native with Python bindings.

   Python users may find Dask more comfortable,
   but Dask is only useful for Python users,
   while Spark can also be used from JVM languages.

-  Dask is one component in the broader Python ecosystem alongside libraries
   like Numpy, Pandas, and Scikit-Learn,
   while Spark is an all-in-one system that re-invents much of the Python world
   in a single package.

   This means that it's often easier to compose Dask with new problem domains,
   but also that you need to install multiple things (like Dask and Pandas or
   Dask and Numpy) rather than just having everything in an all-in-one solution.

-  Apache Spark focuses strongly on traditional business intelligence workloads,
   like ETL, SQL queries, and then some lightweight machine learning,
   while Dask is more general purpose.

   This means that Dask is much more flexible and can handle other problem
   domains like multi-dimensional arrays, GIS, advanced machine learning, and
   custom systems, but that it is less focused and less tuned on typical SQL
   style computations.

   If you mostly want to focus on SQL queries then Spark is probably a better
   bet.  If you want to support a wide variety of custom workloads then Dask
   might be more natural.

See the section :doc:`spark`.


Are there companies that we can get support from?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several companies that offer support for dask in different capacities. See
`Paid support <https://docs.dask.org/en/latest/support.html#paid-support>`_ for a full list.


For IT
------


How would I set up Dask on institutional hardware?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You already have cluster resources.
Dask can run on them today without significant change.

Most institutional clusters today have a resource manager.
This is typically managed by IT, with some mild permissions given to users to
launch jobs.  Dask works with all major resource managers today, including
those on Hadoop, HPC, Kubernetes, and Cloud clusters.

1.  **Hadoop/Spark**: If you have a Hadoop/Spark cluster, such as one purchased
    through Cloudera/Hortonworks/MapR then you will likely want to deploy Dask
    with YARN, the resource manager that deploys services like Hadoop, Spark,
    Hive, and others.

    To help with this, you'll likely want to use `Dask-Yarn <https://yarn.dask.org>`_.

2.  **HPC**: If you have an HPC machine that runs resource managers like SGE,
    SLURM, PBS, LSF, Torque, Condor, or other job batch queuing systems, then
    users can launch Dask on these systems today using either:

    - `Dask Jobqueue <https://jobqueue.dask.org>`_ , which uses typical
      ``qsub``, ``sbatch``, ``bsub`` or other submission tools in interactive
      settings.
    - `Dask MPI <https://mpi.dask.org>`_ which uses MPI for deployment in
      batch settings

    For more information see :doc:`how-to/deploy-dask/hpc`

3.  **Kubernetes/Cloud**: Newer clusters may employ Kubernetes for deployment.
    This is particularly commonly used today on major cloud providers,
    all of which provide hosted Kubernetes as a service.  People today use Dask
    on Kubernetes using either of the following:

    - **Helm**: an easy way to stand up a long-running Dask cluster and
      Jupyter notebook

    - **Dask-Kubernetes**: for native Kubernetes integration for fast moving
      or ephemeral deployments.

    For more information see :doc:`how-to/deploy-dask/kubernetes`


Is Dask secure?
~~~~~~~~~~~~~~~

Dask is deployed today within highly secure institutions,
including major financial, healthcare, and government agencies.

That being said it's worth noting that, by it's very nature, Dask enables the
execution of arbitrary user code on a large set of machines. Care should be
taken to isolate, authenticate, and govern access to these machines.  Fortunately,
your institution likely already does this and uses standard technologies like
SSL/TLS, Kerberos, and other systems with which Dask can integrate.


Do I need to purchase a new cluster?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No.  It is easy to run Dask today on most clusters.
If you have a pre-existing HPC or Spark/Hadoop cluster then that will be fine
to start running Dask.

You can start using Dask without any capital expenditure.


How do I manage users?
~~~~~~~~~~~~~~~~~~~~~~

Dask doesn't manage users, you likely have existing systems that do this well.
In a large institutional setting we assume that you already have a resource
manager like Yarn (Hadoop), Kubernetes, or PBS/SLURM/SGE/LSF/..., each of which
have excellent user management capabilities, which are likely preferred by your
IT department anyway.

Dask is designed to operate with user-level permissions, which means that
your data science users should be able to ask those systems mentioned above for
resources, and have their processes tracked accordingly.

However, there are institutions where analyst-level users aren't given direct access to
the cluster.  This is particularly common in Cloudera/Hortonworks Hadoop/Spark deployments.
In these cases some level of explicit indirection may be required.  For this, we
recommend the `Dask Gateway project <https://gateway.dask.org>`_, which uses IT-level
permissions to properly route authenticated users into secure resources.


How do I manage software environments?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This depends on your cluster resource manager:

-  Most HPC users use their network file system
-  Hadoop/Spark/Yarn users package their environment into a tarball and ship it
   around with HDFS (Dask-Yarn integrates with `Conda Pack
   <https://conda.github.io/conda-pack/>`_ for this capability)
-  Kubernetes or Cloud users use Docker images

In each case Dask integrates with existing processes and technologies
that are well understood and familiar to the institution.


How does Dask communicate data between machines?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask usually communicates over TCP, using msgpack for small administrative
messages, and its own protocol for efficiently passing around large data.
The scheduler and each worker host their own TCP server, making Dask a
distributed peer-to-peer network that uses point-to-point communication.
We do not use Spark-style shuffle systems.  We do not use MPI-style
collectives.  Everything is direct point-to-point.

For high performance networks you can use either TCP-over-Infiniband for about
1 GB/s bandwidth, or UCX (experimental) for full speed communication.


Are deployments long running, or ephemeral?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We see both, but ephemeral deployments are more common.

Most Dask use today is about enabling data science or data engineering users to
scale their interactive workloads across the cluster.
These are typically either interactive sessions with Jupyter, or batch scripts
that run at a pre-defined time.  In both cases, the user asks the resource
manager for a bunch of machines, does some work, and then gives up those
machines.

Some institutions also use Dask in an always-on fashion, either handling
real-time traffic in a scalable way, or responding to a broad set of
interactive users with large datasets that it keeps resident in memory.


For Users
---------

Will Dask "just work" on our existing code?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No, you will need to make modifications,
but these modifications are usually small.

The vast majority of lines of business logic within your institution
will not have to change, assuming that they are in Python and use tooling like
Numpy, Pandas and Scikit-Learn.

How well does Dask scale?  What are Dask's limitations?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The largest Dask deployments that we see today are on around 1000 multi-core
machines, perhaps 20,000 cores in total, but these are rare.
Most institutional-level problems (1-100 TB) are well solved by deployments of 10-50 nodes.

Technically, the back-of-the-envelope number to keep in mind is that each task
(an individual Python function call) in Dask has an overhead of around *200
microseconds*.  So if these tasks take 1 second each, then Dask can saturate
around 5000 cores before scheduling overhead dominates costs.  As workloads
reach this limit they are encouraged to use larger chunk sizes to compensate.
The *vast majority* of institutional users though do not reach this limit.
For more information you may want to peruse our :doc:`best practices
<best-practices>`

Is Dask resilient?  What happens when a machine goes down?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes, Dask is resilient to the failure of worker nodes.  It knows how it came to
any result, and can replay the necessary work on other machines if one goes
down.

If Dask's centralized scheduler goes down then you would need to resubmit the
computation.  This is a fairly standard level of resiliency today, shared with
other tooling like Apache Spark, Flink, and others.

The resource managers that host Dask, like Yarn or Kubernetes, typically
provide long-term 24/7 resilience for always-on operation.

Is the API exactly the same as NumPy/Pandas/Scikit-Learn?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No, but it's very close.  That being said your data scientists will still
have to learn some things.

What we find is that the Numpy/Pandas/Scikit-Learn APIs aren't the challenge
when institutions adopt Dask.  When API inconsistencies do exist, even
modestly skilled programmers are able to understand why and work around them
without much pain.

Instead, the challenge is building intuition around parallel performance.
We've all built up a mental model for what is fast and slow on a single
machine.  This model changes when we factor in network communication and
parallel algorithms, and the performance that we get for familiar operations
can be surprising.

Our main solution to build this intuition, other than
accumulated experience, is Dask's :doc:`Diagnostic Dashboard
<diagnostics-distributed>`.
The dashboard delivers a ton of visual feedback to users as they are running
their computation to help them understand what is going on.  This both helps
them to identify and resolve immediate bottlenecks, and also builds up that
parallel performance intuition surprisingly quickly.


How much performance tuning does Dask require?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*Some other systems are notoriously hard to tune for optimal performance.
What is Dask's story here?  How many knobs are there that we need to be aware
of?*

Like the rest of the Python software tools, Dask puts a lot of effort into
having sane defaults.  Dask workers automatically detect available memory and
cores, and choose sensible defaults that are decent in most situations.  Dask
algorithms similarly provide decent choices by default, and informative warnings
when tricky situations arise, so that, in common cases, things should be fine.

The most common knobs to tune include the following:

-   The thread/process mixture to deal with GIL-holding computations (which are
    rare in Numpy/Pandas/Scikit-Learn workflows)
-   Partition size, like if should you have 100 MB chunks or 1 GB chunks

That being said, almost no institution's needs are met entirely by the common
case, and given the variety of problems that people throw at Dask,
exceptional problems are commonplace.
In these cases we recommend watching the dashboard during execution to see what
is going on.  It can commonly inform you what's going wrong, so that you can
make changes to your system.


What Data formats does Dask support?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because Dask builds on NumPy and Pandas, it supports most formats that they
support, which is most formats.
That being said, not all formats are well suited for
parallel access.  In general people using the following formats are usually
pretty happy:

-  **Tabular:** Parquet, ORC, CSV, Line Delimited JSON, Avro, text
-  **Arrays:** HDF5, NetCDF, Zarr, GRIB

More generally, if you have a Python function that turns a chunk of your stored
data into a Pandas dataframe or Numpy array then Dask can probably call that
function many times without much effort.

For groups looking for advice on which formats to use, we recommend Parquet
for tables and Zarr or HDF5 for arrays.


Does Dask have a SQL interface?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask supports various ways to communicate with SQL databases, some
requiring extra packages to be installed; see the section
:doc:`dataframe-sql`.


Does Dask work on GPUs?
~~~~~~~~~~~~~~~~~~~~~~~

Yes! Dask works with GPUs in a few ways.

The `RAPIDS <https://rapids.ai>`_ libraries provide a GPU-accelerated
Pandas-like library,
`cuDF <https://github.com/rapidsai/cudf>`_,
which interoperates well and is tested against Dask DataFrame.

`Chainer's CuPy <https://cupy.chainer.org/>`_ library provides a GPU
accelerated NumPy-like library that interoperates nicely with Dask Array.

For custom workflows people use Dask alongside GPU-accelerated libraries like PyTorch and
TensorFlow to manage workloads across several machines.  They typically use
Dask's custom APIs, notably :doc:`Delayed <delayed>` and :doc:`Futures
<futures>`.

See the section :doc:`gpu`.


How do I cite Dask?
~~~~~~~~~~~~~~~~~~~

Dask is developed by many people from many institutions.  Some of these
developers are academics who depend on academic citations to justify their
efforts.  Unfortunately, no single citation can do all of these developers (and
the developers to come) sufficient justice.  Instead, we choose to use a single
blanket citation for all developers past and present.

To cite Dask in publications, please use the following::

   Dask Development Team (2016). Dask: Library for dynamic task scheduling
   URL https://dask.org

A BibTeX entry for LaTeX users follows::

   @Manual{,
     title = {Dask: Library for dynamic task scheduling},
     author = {{Dask Development Team}},
     year = {2016},
     url = {https://dask.org},
   }

The full author list is available using ``git`` (e.g. ``git shortlog -ns``).


Are there papers about Dask?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rocklin, Matthew. "Dask: Parallel Computation with Blocked algorithms and Task
Scheduling." (2015).
`PDF <https://conference.scipy.org/proceedings/scipy2015/pdfs/matthew_rocklin.pdf>`_.

::

   @InProceedings{ matthew_rocklin-proc-scipy-2015,
     author    = { Matthew Rocklin },
     title     = { Dask: Parallel Computation with Blocked algorithms and Task Scheduling },
     booktitle = { Proceedings of the 14th Python in Science Conference },
     pages     = { 130 - 136 },
     year      = { 2015 },
     editor    = { Kathryn Huff and James Bergstra }
   }


For Marketing
-------------

There is a special subsite dedicated to addressing marketing concerns. You can
find it at `marketing.dask.org <https://marketing.dask.org>`_.

Do you have any standardized logos?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes! You can find them at :doc:`logos`.
:orphan:

Diagnostics (local)
====================

Profiling parallel code can be challenging, but ``dask.diagnostics`` provides
functionality to aid in profiling and inspecting execution with the
:doc:`local task scheduler <scheduling>`.

This page describes the following few built-in options:

1.  ProgressBar
2.  Profiler
3.  ResourceProfiler
4.  CacheProfiler

Furthermore, this page then provides instructions on how to build your own custom diagnostic.

.. currentmodule:: dask.diagnostics


Progress Bar
------------

.. autosummary::
   ProgressBar

The ``ProgressBar`` class builds on the scheduler callbacks described above to
display a progress bar in the terminal or notebook during computation. This can
give a nice feedback during long running graph execution. It can be used as a
context manager around calls to ``get`` or ``compute`` to profile the
computation:

.. code-block:: python

    >>> import dask.array as da
    >>> from dask.diagnostics import ProgressBar
    >>> a = da.random.normal(size=(10000, 10000), chunks=(1000, 1000))
    >>> res = a.dot(a.T).mean(axis=0)

    >>> with ProgressBar():
    ...     out = res.compute()
    [########################################] | 100% Completed | 17.1 s

or registered globally using the ``register`` method:

.. code-block:: python

    >>> pbar = ProgressBar()
    >>> pbar.register()
    >>> out = res.compute()
    [########################################] | 100% Completed | 17.1 s

To unregister from the global callbacks, call the ``unregister`` method:

.. code-block:: python

    >>> pbar.unregister()



Profiler
--------

.. autosummary::
   Profiler

Dask provides a few tools for profiling execution. As with the ``ProgressBar``,
they each can be used as context managers or registered globally.

The ``Profiler`` class is used to profile Dask's execution at the task level.
During execution, it records the following information for each task:

1. Key
2. Task
3. Start time in seconds since the epoch
4. Finish time in seconds since the epoch
5. Worker id


ResourceProfiler
----------------

.. autosummary::
   ResourceProfiler

The ``ResourceProfiler`` class is used to profile Dask's execution at the
resource level. During execution, it records the following information
for each timestep:

1. Time in seconds since the epoch
2. Memory usage in MB
3. % CPU usage

The default timestep is 1 second, but can be set manually using the ``dt``
keyword:

.. code-block:: python

    >>> from dask.diagnostics import ResourceProfiler
    >>> rprof = ResourceProfiler(dt=0.5)


CacheProfiler
-------------

.. autosummary::
   CacheProfiler

The ``CacheProfiler`` class is used to profile Dask's execution at the scheduler
cache level. During execution, it records the following information for each
task:

1. Key
2. Task
3. Size metric
4. Cache entry time in seconds since the epoch
5. Cache exit time in seconds since the epoch

Here the size metric is the output of a function called on the result of each
task. The default metric is to count each task (``metric`` is 1 for all tasks).
Other functions may be used as a metric instead through the ``metric`` keyword.
For example, the ``nbytes`` function found in ``cachey`` can be used to measure
the number of bytes in the scheduler cache:

.. code-block:: python

    >>> from dask.diagnostics import CacheProfiler
    >>> from cachey import nbytes
    >>> cprof = CacheProfiler(metric=nbytes)


Example
-------

As an example to demonstrate using the diagnostics, we'll profile some linear
algebra done with Dask Array. We'll create a random array, take its QR
decomposition, and then reconstruct the initial array by multiplying the Q and
R components together. Note that since the profilers (and all diagnostics) are
just context managers, multiple profilers can be used in a with block:

.. code-block:: python

    >>> import dask.array as da
    >>> from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler
    >>> a = da.random.random(size=(10000, 1000), chunks=(1000, 1000))
    >>> q, r = da.linalg.qr(a)
    >>> a2 = q.dot(r)

    >>> with Profiler() as prof, ResourceProfiler(dt=0.25) as rprof,
    ...         CacheProfiler() as cprof:
    ...     out = a2.compute()

The results of each profiler are stored in their ``results`` attribute as a
list of ``namedtuple`` objects:

.. code-block:: python

    >>> prof.results[0]
    TaskData(key=('tsqr-8d16e396b237bf7a731333130d310cb9_QR_st1', 5, 0),
             task=(qr, (_apply_random, 'random_sample', 1060164455, (1000, 1000), (), {})),
             start_time=1454368444.493292,
             end_time=1454368444.902987,
             worker_id=4466937856)

    >>> rprof.results[0]
    ResourceData(time=1454368444.078748, mem=74.100736, cpu=0.0)

    >>> cprof.results[0]
    CacheData(key=('tsqr-8d16e396b237bf7a731333130d310cb9_QR_st1', 7, 0),
              task=(qr, (_apply_random, 'random_sample', 1310656009, (1000, 1000), (), {})),
              metric=1,
              cache_time=1454368444.49662,
              free_time=1454368446.769452)

These can be analyzed separately or viewed in a bokeh plot using the provided
``visualize`` method on each profiler:

.. code-block:: python

    >>> prof.visualize()


.. raw:: html

    <iframe src="_static/profile.html"
            marginwidth="0" marginheight="0" scrolling="no"
            width="650" height="300" style="border:none"></iframe>

To view multiple profilers at the same time, the :func:`dask.diagnostics.visualize`
function can be used. This takes a list of profilers and creates a vertical
stack of plots aligned along the x-axis:

.. code-block:: python

    >>> from dask.diagnostics import visualize
    >>> visualize([prof, rprof, cprof])


.. raw:: html

    <iframe src="_static/stacked_profile.html"
            marginwidth="0" marginheight="0" scrolling="no"
            width="650" height="700" style="border:none"></iframe>


Looking at the above figure, from top to bottom:

1. The results from the ``Profiler`` object: This shows the execution time for
   each task as a rectangle, organized along the y-axis by worker (in this case
   threads). Similar tasks are grouped by color and, by hovering over each task,
   one can see the key and task that each block represents.

2. The results from the ``ResourceProfiler`` object: This shows two lines, one
   for total CPU percentage used by all the workers, and one for total memory
   usage.

3. The results from the ``CacheProfiler`` object: This shows a line for each
   task group, plotting the sum of the current ``metric`` in the cache against
   time. In this case it's the default metric (count) and the lines represent
   the number of each object in the cache at time. Note that the grouping and
   coloring is the same as for the ``Profiler`` plot, and that the task
   represented by each line can be found by hovering over the line.

From these plots we can see that the initial tasks (calls to
``numpy.random.random`` and ``numpy.linalg.qr`` for each chunk) are run
concurrently, but only use slightly more than 100\% CPU. This is because the
call to ``numpy.linalg.qr`` currently doesn't release the Global Interpreter
Lock (GIL), so those calls can't truly be done in parallel. Next, there's a reduction
step where all the blocks are combined. This requires all the results from the
first step to be held in memory, as shown by the increased number of results in
the cache, and increase in memory usage. Immediately after this task ends, the
number of elements in the cache decreases, showing that they were only needed
for this step. Finally, there's an interleaved set of calls to ``dot`` and
``sum``. Looking at the CPU plot, it shows that these run both concurrently and in
parallel, as the CPU percentage spikes up to around 350\%.


Custom Callbacks
----------------

.. autosummary:: Callback

Schedulers based on ``dask.local.get_async`` (currently
``dask.get``, ``dask.threaded.get``, and ``dask.multiprocessing.get``)
accept five callbacks, allowing for inspection of scheduler execution.

The callbacks are:

1. ``start(dsk)``: Run at the beginning of execution, right before the
state is initialized.  Receives the Dask graph

2. ``start_state(dsk, state)``: Run at the beginning of execution, right
after the state is initialized.  Receives the Dask graph and scheduler state

3. ``pretask(key, dsk, state)``: Run every time a new task is started.
Receives the key of the task to be run, the Dask graph, and the scheduler state

4. ``posttask(key, result, dsk, state, id)``: Run every time a task is finished.
Receives the key of the task that just completed, the result, the Dask graph,
the scheduler state, and the id of the worker that ran the task

5. ``finish(dsk, state, errored)``: Run at the end of execution, right before the
result is returned. Receives the Dask graph, the scheduler state, and a boolean
indicating whether or not the exit was due to an error

Custom diagnostics can be created either by instantiating the ``Callback``
class with the some of the above methods as keywords or by subclassing the
``Callback`` class.
Here we create a class that prints the name of every key as it's computed:

.. code-block:: python

    from dask.callbacks import Callback
    class PrintKeys(Callback):
        def _pretask(self, key, dask, state):
            """Print the key of every task as it's started"""
            print("Computing: {0}!".format(repr(key)))

This can now be used as a context manager during computation:

.. code-block:: python

    >>> from operator import add, mul
    >>> dsk = {'a': (add, 1, 2), 'b': (add, 3, 'a'), 'c': (mul, 'a', 'b')}

    >>> with PrintKeys():
    ...     get(dsk, 'c')
    Computing 'a'!
    Computing 'b'!
    Computing 'c'!

Alternatively, functions may be passed in as keyword arguments to ``Callback``:

.. code-block:: python

    >>> def printkeys(key, dask, state):
    ...    print("Computing: {0}!".format(repr(key)))

    >>> with Callback(pretask=printkeys):
    ...     get(dsk, 'c')
    Computing 'a'!
    Computing 'b'!
    Computing 'c'!


API
---

.. autosummary::
   CacheProfiler
   Callback
   Profiler
   ProgressBar
   ResourceProfiler
   visualize

.. autofunction:: ProgressBar
.. autofunction:: Profiler
.. autofunction:: ResourceProfiler
.. autofunction:: CacheProfiler
.. autofunction:: Callback
.. autofunction:: visualize
.. _graph_manipulation:

Advanced graph manipulation
===========================
There are some situations where computations with Dask collections will result in
suboptimal memory usage (e.g. an entire Dask DataFrame is loaded into memory).
This may happen when Dask’s scheduler doesn’t automatically delay the computation of
nodes in a task graph to avoid occupying memory with their output for prolonged periods
of time, or in scenarios where recalculating nodes is much cheaper than holding their
output in memory.

This page highlights a set of graph manipulation utilities which can be used to help
avoid these scenarios. In particular, the utilities described below rewrite the
underlying Dask graph for Dask collections, producing equivalent collections with
different sets of keys.

Consider the following example:

.. code-block:: python

   >>> import dask.array as da
   >>> x = da.random.normal(size=500_000_000, chunks=100_000)
   >>> x_mean = x.mean()
   >>> y = (x - x_mean).max().compute()

The above example computes the largest value of a distribution after removing its bias.
This involves loading the chunks of ``x`` into memory in order to compute ``x_mean``.
However, since the ``x`` array is needed later in the computation to compute ``y``, the
entire ``x`` array is kept in memory. For large Dask Arrays this can be very
problematic.

To alleviate the need for the entire ``x`` array to be kept in memory, one could rewrite
the last line as follows:

.. code-block:: python

   >>> from dask.graph_manipulation import bind
   >>> xb = bind(x, x_mean)
   >>> y = (xb - x_mean).max().compute()

Here we use :func:`~dask.graph_manipulation.bind` to create a new Dask Array, ``xb``,
which produces exactly the same output as ``x``, but whose underlying Dask graph has
different keys than ``x``, and will only be computed after ``x_mean`` has been
calculated.

This results in the chunks of ``x`` being computed and immediately individually reduced
by ``mean``; then recomputed and again immediately pipelined into the subtraction
followed by reduction with ``max``. This results in a much smaller peak memory usage as
the full ``x`` array is no longer loaded into memory. However, the tradeoff is that the
compute time increases as ``x`` is computed twice.


API
---

.. currentmodule:: dask.graph_manipulation

.. autosummary::

   checkpoint
   wait_on
   bind
   clone


Definitions
~~~~~~~~~~~

.. autofunction:: checkpoint
.. autofunction:: wait_on
.. autofunction:: bind
.. autofunction:: clone
:orphan:

Visualize task graphs
---------------------

.. currentmodule:: dask

.. autosummary::
   visualize

Before executing your computation you might consider visualizing the underlying task graph.
By looking at the inter-connectedness of tasks
you can learn more about potential bottlenecks
where parallelism may not be possible,
or areas where many tasks depend on each other,
which may cause a great deal of communication.

The ``.visualize`` method and ``dask.visualize`` function work exactly like
the ``.compute`` method and ``dask.compute`` function,
except that rather than computing the result,
they produce an image of the task graph.

By default the task graph is rendered from top to bottom.
In the case that you prefer to visualize it from left to right, pass
``rankdir="LR"`` as a keyword argument to ``.visualize``.

.. code-block:: python

   import dask.array as da
   x = da.ones((15, 15), chunks=(5, 5))

   y = x + x.T

   # y.compute()
   y.visualize(filename='transpose.svg')

.. image:: images/transpose.svg
   :alt: Dask task graph for adding an array to its transpose

Note that the ``visualize`` function is powered by the `GraphViz <https://www.graphviz.org/>`_
system library.  This library has a few considerations:

1.  You must install both the graphviz system library (with tools like apt-get, yum, or brew)
    *and* the graphviz Python library.
    If you use Conda then you need to install ``python-graphviz``,
    which will bring along the ``graphviz`` system library as a dependency.
2.  Graphviz takes a while on graphs larger than about 100 nodes.
    For large computations you might have to simplify your computation a bit
    for the visualize method to work well.

Dask Internals
==============

This section is intended for contributors and power users who are interested in
learning more about how Dask works internally.

.. toctree::
   :maxdepth: 1

   user-interfaces.rst
   understanding-performance.rst
   phases-of-computation.rst
   order.rst
   caching.rst
   shared.rst
   scheduling-policy.rst
Array
=====

.. toctree::
   :maxdepth: 1
   :hidden:

   array-best-practices.rst
   array-chunks.rst
   array-creation.rst
   array-overlap.rst
   array-design.rst
   array-sparse.rst
   array-stats.rst
   array-linear-operator.rst
   array-slicing.rst
   array-assignment.rst
   array-stack.rst
   array-gufunc.rst

Dask Array implements a subset of the NumPy ndarray interface using blocked
algorithms, cutting up the large array into many small arrays. This lets us
compute on arrays larger than memory using all of our cores.  We coordinate
these blocked algorithms using Dask graphs.

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/9h_61hXCDuI"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

Examples
--------

Visit https://examples.dask.org/array.html to see and run examples using
Dask Array.

Design
------

.. image:: images/dask-array.svg
   :alt: Dask arrays coordinate many numpy arrays
   :align: right
   :scale: 35%

Dask arrays coordinate many NumPy arrays (or "duck arrays" that are
sufficiently NumPy-like in API such as CuPy or Sparse arrays) arranged into a
grid. These arrays may live on disk or on other machines.

New duck array chunk types (types below Dask on
`NEP-13's type-casting hierarchy`_) can be registered via
:func:`~dask.array.register_chunk_type`. Any other duck array types that are
not registered will be deferred to in binary operations and NumPy
ufuncs/functions (that is, Dask will return ``NotImplemented``). Note, however,
that *any* ndarray-like type can be inserted into a Dask Array using
:func:`~dask.array.Array.from_array`.

Common Uses
-----------

Dask Array is used in fields like atmospheric and oceanographic science, large
scale imaging, genomics, numerical algorithms for optimization or statistics,
and more.

Scope
-----

Dask arrays support most of the NumPy interface like the following:

-  Arithmetic and scalar mathematics: ``+, *, exp, log, ...``
-  Reductions along axes: ``sum(), mean(), std(), sum(axis=0), ...``
-  Tensor contractions / dot products / matrix multiply: ``tensordot``
-  Axis reordering / transpose: ``transpose``
-  Slicing: ``x[:100, 500:100:-2]``
-  Fancy indexing along single axes with lists or NumPy arrays: ``x[:, [10, 1, 5]]``
-  Array protocols like ``__array__`` and ``__array_ufunc__``
-  Some linear algebra: ``svd, qr, solve, solve_triangular, lstsq``
-  ...


However, Dask Array does not implement the entire NumPy interface.  Users expecting this
will be disappointed.  Notably, Dask Array lacks the following features:

-   Much of ``np.linalg`` has not been implemented.
    This has been done by a number of excellent BLAS/LAPACK implementations,
    and is the focus of numerous ongoing academic research projects
-   Arrays with unknown shapes do not support all operations
-   Operations like ``sort`` which are notoriously
    difficult to do in parallel, and are of somewhat diminished value on very
    large data (you rarely actually need a full sort).
    Often we include parallel-friendly alternatives like ``topk``
-   Dask Array doesn't implement operations like ``tolist`` that would be very
    inefficient for larger datasets. Likewise, it is very inefficient to iterate
    over a Dask array with for loops
-   Dask development is driven by immediate need, hence many lesser used
    functions have not been implemented. Community contributions are encouraged

See :doc:`the dask.array API<array-api>` for a more extensive list of
functionality.

Execution
---------

By default, Dask Array uses the threaded scheduler in order to avoid data
transfer costs, and because NumPy releases the GIL well.  It is also quite
effective on a cluster using the `dask.distributed`_ scheduler.

.. _`dask.distributed`: https://distributed.dask.org/en/latest/
.. _`NEP-13's type-casting hierarchy`: https://numpy.org/neps/nep-0013-ufunc-overrides.html#type-casting-hierarchy
====
Dask
====

*Dask is a flexible library for parallel computing in Python.*

Dask is composed of two parts:

1.  **Dynamic task scheduling** optimized for computation. This is similar to
    *Airflow, Luigi, Celery, or Make*, but optimized for interactive
    computational workloads.
2.  **"Big Data" collections** like parallel arrays, dataframes, and lists that
    extend common interfaces like *NumPy, Pandas, or Python iterators* to
    larger-than-memory or distributed environments. These parallel collections
    run on top of dynamic task schedulers.

Dask emphasizes the following virtues:

*  **Familiar**: Provides parallelized NumPy array and Pandas DataFrame objects
*  **Flexible**: Provides a task scheduling interface for more custom workloads
   and integration with other projects.
*  **Native**: Enables distributed computing in pure Python with access to
   the PyData stack.
*  **Fast**: Operates with low overhead, low latency, and minimal serialization
   necessary for fast numerical algorithms
*  **Scales up**: Runs resiliently on clusters with 1000s of cores
*  **Scales down**: Trivial to set up and run on a laptop in a single process
*  **Responsive**: Designed with interactive computing in mind, it provides rapid
   feedback and diagnostics to aid humans

|

.. figure:: images/dask-overview.svg
   :alt: Dask is composed of three parts. "Collections" create "Task Graphs" which are then sent to the "Scheduler" for execution. There are two types of schedulers that are described in more detail below.
   :align: center

   High level collections are used to generate task graphs which can be executed by schedulers on a single machine or a cluster.

|

See the `dask.distributed documentation (separate website)
<https://distributed.dask.org/en/latest/>`_ for more technical information
on Dask's distributed scheduler.

Familiar user interface
-----------------------

**Dask DataFrame** mimics Pandas - :doc:`documentation <dataframe>`

.. code-block:: python

    import pandas as pd                     import dask.dataframe as dd
    df = pd.read_csv('2015-01-01.csv')      df = dd.read_csv('2015-*-*.csv')
    df.groupby(df.user_id).value.mean()     df.groupby(df.user_id).value.mean().compute()

**Dask Array** mimics NumPy - :doc:`documentation <array>`

.. code-block:: python

   import numpy as np                       import dask.array as da
   f = h5py.File('myfile.hdf5')             f = h5py.File('myfile.hdf5')
   x = np.array(f['/small-data'])           x = da.from_array(f['/big-data'],
                                                              chunks=(1000, 1000))
   x - x.mean(axis=1)                       x - x.mean(axis=1).compute()

**Dask Bag** mimics iterators, Toolz, and PySpark - :doc:`documentation <bag>`

.. code-block:: python

   import dask.bag as db
   b = db.read_text('2015-*-*.json.gz').map(json.loads)
   b.pluck('name').frequencies().topk(10, lambda pair: pair[1]).compute()

**Dask Delayed** mimics for loops and wraps custom code - :doc:`documentation <delayed>`

.. code-block:: python

   from dask import delayed
   L = []
   for fn in filenames:                  # Use for loops to build up computation
       data = delayed(load)(fn)          # Delay execution of function
       L.append(delayed(process)(data))  # Build connections between variables

   result = delayed(summarize)(L)
   result.compute()

The **concurrent.futures** interface provides general submission of custom
tasks: - :doc:`documentation <futures>`

.. code-block:: python

   from dask.distributed import Client
   client = Client('scheduler:port')

   futures = []
   for fn in filenames:
       future = client.submit(load, fn)
       futures.append(future)

   summary = client.submit(summarize, futures)
   summary.result()


Scales from laptops to clusters
-------------------------------

Dask is convenient on a laptop.  It :doc:`installs <install>` trivially with
``conda`` or ``pip`` and extends the size of convenient datasets from "fits in
memory" to "fits on disk".

Dask can scale to a cluster of 100s of machines. It is resilient, elastic, data
local, and low latency.  For more information, see the documentation about the
`distributed scheduler`_.

This ease of transition between single-machine to moderate cluster enables
users to both start simple and grow when necessary.


Complex Algorithms
------------------

Dask represents parallel computations with :doc:`task graphs<graphs>`. These
directed acyclic graphs may have arbitrary structure, which enables both
developers and users the freedom to build sophisticated algorithms and to
handle messy situations not easily managed by the ``map/filter/groupby``
paradigm common in most data engineering frameworks.

We originally needed this complexity to build complex algorithms for
n-dimensional arrays but have found it to be equally valuable when dealing with
messy situations in everyday problems.

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Getting Started

   install.rst
   10-minutes-to-dask.rst
   presentations.rst
   best-practices.rst
   how-to/index.rst
   faq.rst

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Fundamentals

   array.rst
   bag.rst
   dataframe.rst
   delayed.rst
   futures.rst
   scheduling.rst
   graphs.rst
   deploying.rst
   internals.rst

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Reference

   api.rst
   develop.rst
   changelog.rst
   configuration.rst

.. _`Anaconda Inc`: https://www.anaconda.com
.. _`3-clause BSD license`: https://github.com/dask/dask/blob/main/LICENSE.txt

.. _`#dask tag`: https://stackoverflow.com/questions/tagged/dask
.. _`GitHub issue tracker`: https://github.com/dask/dask/issues
.. _`xarray`: https://xarray.pydata.org/en/stable/
.. _`scikit-image`: https://scikit-image.org/docs/stable/
.. _`scikit-allel`: https://scikits.appspot.com/scikit-allel
.. _`pandas`: https://pandas.pydata.org/pandas-docs/version/0.17.0/
.. _`distributed scheduler`: https://distributed.dask.org/en/latest/
Sparse Arrays
=============

By swapping out in-memory NumPy arrays with in-memory sparse arrays, we can
reuse the blocked algorithms of Dask's Array to achieve parallel and distributed
sparse arrays.

The blocked algorithms in Dask Array normally parallelize around in-memory
NumPy arrays.  However, if another in-memory array library supports the NumPy
interface, then it too can take advantage of Dask Array's parallel algorithms.
In particular the `sparse <https://github.com/pydata/sparse/>`_ array library
satisfies a subset of the NumPy API and works well with (and is tested against)
Dask Array.

Example
-------

Say we have a Dask array with mostly zeros:

.. code-block:: python

   x = da.random.random((100000, 100000), chunks=(1000, 1000))
   x[x < 0.95] = 0

We can convert each of these chunks of NumPy arrays into a **sparse.COO** array:

.. code-block:: python

   import sparse
   s = x.map_blocks(sparse.COO)

Now, our array is not composed of many NumPy arrays, but rather of many
sparse arrays.  Semantically, this does not change anything.  Operations that
work will continue to work identically (assuming that the behavior of ``numpy`` 
and ``sparse`` are identical), but performance characteristics and storage costs 
may change significantly:

.. code-block:: python

   >>> s.sum(axis=0)[:100].compute()
   <COO: shape=(100,), dtype=float64, nnz=100>

   >>> _.todense()
   array([ 4803.06859272,  4913.94964525,  4877.13266438,  4860.7470773 ,
           4938.94446802,  4849.51326473,  4858.83977856,  4847.81468485,
           ... ])

Requirements
------------

Any in-memory library that copies the NumPy ndarray interface should work here.
The `sparse <https://github.com/pydata/sparse/>`_ library is a minimal
example.  In particular, an in-memory library should implement at least the
following operations:

1.  Simple slicing with slices, lists, and elements (for slicing, rechunking,
    reshaping, etc)
2.  A ``concatenate`` function matching the interface of ``np.concatenate``.
    This must be registered in ``dask.array.core.concatenate_lookup``
3.  All ufuncs must support the full ufunc interface, including ``dtype=`` and
    ``out=`` parameters (even if they don't function properly)
4.  All reductions must support the full ``axis=`` and ``keepdims=`` keywords
    and behave like NumPy in this respect
5.  The array class should follow the ``__array_priority__`` protocol and be
    prepared to respond to other arrays of lower priority
6.  If ``dot`` support is desired, a ``tensordot`` function matching the
    interface of ``np.tensordot`` should be registered in
    ``dask.array.core.tensordot_lookup``

The implementation of other operations like reshape, transpose, etc.,
should follow standard NumPy conventions regarding shape and dtype.  Not
implementing these is fine; the parallel ``dask.array`` will err at runtime if
these operations are attempted.


Mixed Arrays
------------

Dask's Array supports mixing different kinds of in-memory arrays.  This relies
on the in-memory arrays knowing how to interact with each other when necessary.
When two arrays interact, the functions from the array with the highest
``__array_priority__`` will take precedence (for example, for concatenate,
tensordot, etc.).
API
---

.. currentmodule:: dask.dataframe

Dataframe
~~~~~~~~~

.. autosummary::
    :toctree: generated/

    DataFrame
    DataFrame.abs
    DataFrame.add
    DataFrame.align
    DataFrame.all
    DataFrame.any
    DataFrame.append
    DataFrame.apply
    DataFrame.applymap
    DataFrame.assign
    DataFrame.astype
    DataFrame.bfill
    DataFrame.categorize
    DataFrame.columns
    DataFrame.compute
    DataFrame.copy
    DataFrame.corr
    DataFrame.count
    DataFrame.cov
    DataFrame.cummax
    DataFrame.cummin
    DataFrame.cumprod
    DataFrame.cumsum
    DataFrame.describe
    DataFrame.diff
    DataFrame.div
    DataFrame.divide
    DataFrame.drop
    DataFrame.drop_duplicates
    DataFrame.dropna
    DataFrame.dtypes
    DataFrame.eq
    DataFrame.eval
    DataFrame.explode
    DataFrame.ffill
    DataFrame.fillna
    DataFrame.first
    DataFrame.floordiv
    DataFrame.ge
    DataFrame.get_partition
    DataFrame.groupby
    DataFrame.gt
    DataFrame.head
    DataFrame.idxmax
    DataFrame.idxmin
    DataFrame.iloc
    DataFrame.index
    DataFrame.info
    DataFrame.isin
    DataFrame.isna
    DataFrame.isnull
    DataFrame.items
    DataFrame.iterrows
    DataFrame.itertuples
    DataFrame.join
    DataFrame.known_divisions
    DataFrame.last
    DataFrame.le
    DataFrame.loc
    DataFrame.lt
    DataFrame.map_partitions
    DataFrame.mask
    DataFrame.max
    DataFrame.mean
    DataFrame.melt
    DataFrame.memory_usage
    DataFrame.memory_usage_per_partition
    DataFrame.merge
    DataFrame.min
    DataFrame.mod
    DataFrame.mode
    DataFrame.mul
    DataFrame.ndim
    DataFrame.ne
    DataFrame.nlargest
    DataFrame.npartitions
    DataFrame.nsmallest
    DataFrame.partitions
    DataFrame.pivot_table
    DataFrame.pop
    DataFrame.pow
    DataFrame.prod
    DataFrame.quantile
    DataFrame.query
    DataFrame.radd
    DataFrame.random_split
    DataFrame.rdiv
    DataFrame.rename
    DataFrame.repartition
    DataFrame.replace
    DataFrame.resample
    DataFrame.reset_index
    DataFrame.rfloordiv
    DataFrame.rmod
    DataFrame.rmul
    DataFrame.round
    DataFrame.rpow
    DataFrame.rsub
    DataFrame.rtruediv
    DataFrame.sample
    DataFrame.select_dtypes
    DataFrame.sem
    DataFrame.set_index
    DataFrame.shape
    DataFrame.shuffle
    DataFrame.size
    DataFrame.sort_values
    DataFrame.squeeze
    DataFrame.std
    DataFrame.sub
    DataFrame.sum
    DataFrame.tail
    DataFrame.to_bag
    DataFrame.to_csv
    DataFrame.to_dask_array
    DataFrame.to_delayed
    DataFrame.to_hdf
    DataFrame.to_html
    DataFrame.to_json
    DataFrame.to_parquet
    DataFrame.to_records
    DataFrame.to_string
    DataFrame.to_sql
    DataFrame.to_timestamp
    DataFrame.truediv
    DataFrame.values
    DataFrame.var
    DataFrame.visualize
    DataFrame.where

Series
~~~~~~

.. autosummary::
   :toctree: generated/

   Series
   Series.add
   Series.align
   Series.all
   Series.any
   Series.append
   Series.apply
   Series.astype
   Series.autocorr
   Series.between
   Series.bfill
   Series.cat
   Series.clear_divisions
   Series.clip
   Series.clip_lower
   Series.clip_upper
   Series.compute
   Series.copy
   Series.corr
   Series.count
   Series.cov
   Series.cummax
   Series.cummin
   Series.cumprod
   Series.cumsum
   Series.describe
   Series.diff
   Series.div
   Series.drop_duplicates
   Series.dropna
   Series.dt
   Series.dtype
   Series.eq
   Series.explode
   Series.ffill
   Series.fillna
   Series.first
   Series.floordiv
   Series.ge
   Series.get_partition
   Series.groupby
   Series.gt
   Series.head
   Series.idxmax
   Series.idxmin
   Series.isin
   Series.isna
   Series.isnull
   Series.iteritems
   Series.known_divisions
   Series.last
   Series.le
   Series.loc
   Series.lt
   Series.map
   Series.map_overlap
   Series.map_partitions
   Series.mask
   Series.max
   Series.mean
   Series.memory_usage
   Series.memory_usage_per_partition
   Series.min
   Series.mod
   Series.mul
   Series.nbytes
   Series.ndim
   Series.ne
   Series.nlargest
   Series.notnull
   Series.nsmallest
   Series.nunique
   Series.nunique_approx
   Series.persist
   Series.pipe
   Series.pow
   Series.prod
   Series.quantile
   Series.radd
   Series.random_split
   Series.rdiv
   Series.reduction
   Series.repartition
   Series.replace
   Series.rename
   Series.resample
   Series.reset_index
   Series.rolling
   Series.round
   Series.sample
   Series.sem
   Series.shape
   Series.shift
   Series.size
   Series.std
   Series.str
   Series.sub
   Series.sum
   Series.to_bag
   Series.to_csv
   Series.to_dask_array
   Series.to_delayed
   Series.to_frame
   Series.to_hdf
   Series.to_string
   Series.to_timestamp
   Series.truediv
   Series.unique
   Series.value_counts
   Series.values
   Series.var
   Series.visualize
   Series.where


Groupby Operations
~~~~~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe.groupby

DataFrame Groupby
*****************

.. autosummary::
   :toctree: generated/

   DataFrameGroupBy.aggregate
   DataFrameGroupBy.apply
   DataFrameGroupBy.count
   DataFrameGroupBy.cumcount
   DataFrameGroupBy.cumprod
   DataFrameGroupBy.cumsum
   DataFrameGroupBy.get_group
   DataFrameGroupBy.max
   DataFrameGroupBy.mean
   DataFrameGroupBy.min
   DataFrameGroupBy.size
   DataFrameGroupBy.std
   DataFrameGroupBy.sum
   DataFrameGroupBy.var
   DataFrameGroupBy.cov
   DataFrameGroupBy.corr
   DataFrameGroupBy.first
   DataFrameGroupBy.last
   DataFrameGroupBy.idxmin
   DataFrameGroupBy.idxmax
   DataFrameGroupBy.rolling


Series Groupby
**************

.. autosummary::
   :toctree: generated/

   SeriesGroupBy.aggregate
   SeriesGroupBy.apply
   SeriesGroupBy.count
   SeriesGroupBy.cumcount
   SeriesGroupBy.cumprod
   SeriesGroupBy.cumsum
   SeriesGroupBy.get_group
   SeriesGroupBy.max
   SeriesGroupBy.mean
   SeriesGroupBy.min
   SeriesGroupBy.nunique
   SeriesGroupBy.size
   SeriesGroupBy.std
   SeriesGroupBy.sum
   SeriesGroupBy.var
   SeriesGroupBy.first
   SeriesGroupBy.last
   SeriesGroupBy.idxmin
   SeriesGroupBy.idxmax
   SeriesGroupBy.rolling

Custom Aggregation
******************

.. autosummary::
   :toctree: generated/

   Aggregation

Rolling Operations
~~~~~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe

.. autosummary::
   :toctree: generated/

   rolling.map_overlap
   Series.rolling
   DataFrame.rolling

.. currentmodule:: dask.dataframe.rolling

.. autosummary::
   :toctree: generated/

   Rolling.apply
   Rolling.count
   Rolling.kurt
   Rolling.max
   Rolling.mean
   Rolling.median
   Rolling.min
   Rolling.quantile
   Rolling.skew
   Rolling.std
   Rolling.sum
   Rolling.var


Create DataFrames
~~~~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe

.. autosummary::
   :toctree: generated/

   read_csv
   read_table
   read_fwf
   read_parquet
   read_hdf
   read_json
   read_orc
   read_sql_table
   read_sql_query
   read_sql
   from_array
   from_bcolz
   from_dask_array
   from_delayed
   from_pandas

.. currentmodule:: dask.bag

.. autosummary::
   :toctree: generated/

   Bag.to_dataframe

Store DataFrames
~~~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe

.. autosummary::
   :toctree: generated/

   to_csv
   to_parquet
   to_hdf
   to_records
   to_sql
   to_json

Convert DataFrames
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   DataFrame.to_bag
   DataFrame.to_dask_array
   DataFrame.to_delayed

Reshape DataFrames
~~~~~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe.reshape

.. autosummary::
   :toctree: generated/

   get_dummies
   pivot_table
   melt

Concatenate DataFrames
~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe.multi

.. autosummary::
   :toctree: generated/

   DataFrame.merge
   concat
   merge
   merge_asof


Resampling
~~~~~~~~~~

.. currentmodule:: dask.dataframe.tseries.resample

.. autosummary::
   :toctree: generated/

   Resampler
   Resampler.agg
   Resampler.count
   Resampler.first
   Resampler.last
   Resampler.max
   Resampler.mean
   Resampler.median
   Resampler.min
   Resampler.nunique
   Resampler.ohlc
   Resampler.prod
   Resampler.quantile
   Resampler.sem
   Resampler.size
   Resampler.std
   Resampler.sum
   Resampler.var


Dask Metadata
~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe.utils

.. autosummary::
   :toctree: generated/

   make_meta

Other functions
~~~~~~~~~~~~~~~

.. currentmodule:: dask.dataframe

.. autosummary::
   :toctree: generated/

   compute
   map_partitions

   to_datetime
   to_numeric

.. _optimization:

Optimization
============

Performance can be significantly improved in different contexts by making
small optimizations on the Dask graph before calling the scheduler.

The ``dask.optimization`` module contains several functions to transform graphs
in a variety of useful ways. In most cases, users won't need to interact with
these functions directly, as specialized subsets of these transforms are done
automatically in the Dask collections (``dask.array``, ``dask.bag``, and
``dask.dataframe``). However, users working with custom graphs or computations
may find that applying these methods results in substantial speedups.

In general, there are two goals when doing graph optimizations:

1. Simplify computation
2. Improve parallelism

Simplifying computation can be done on a graph level by removing unnecessary
tasks (``cull``), or on a task level by replacing expensive operations with
cheaper ones (``RewriteRule``).

Parallelism can be improved by reducing
inter-task communication, whether by fusing many tasks into one (``fuse``), or
by inlining cheap operations (``inline``, ``inline_functions``).

Below, we show an example walking through the use of some of these to optimize
a task graph.

Example
-------

Suppose you had a custom Dask graph for doing a word counting task:

.. code-block:: python

    >>> def print_and_return(string):
    ...     print(string)
    ...     return string

    >>> def format_str(count, val, nwords):
    ...     return (f'word list has {count} occurrences of '
    ...             f'{val}, out of {nwords} words')

    >>> dsk = {'words': 'apple orange apple pear orange pear pear',
    ...        'nwords': (len, (str.split, 'words')),
    ...        'val1': 'orange',
    ...        'val2': 'apple',
    ...        'val3': 'pear',
    ...        'count1': (str.count, 'words', 'val1'),
    ...        'count2': (str.count, 'words', 'val2'),
    ...        'count3': (str.count, 'words', 'val3'),
    ...        'format1': (format_str, 'count1', 'val1', 'nwords'),
    ...        'format2': (format_str, 'count2', 'val2', 'nwords'),
    ...        'format3': (format_str, 'count3', 'val3', 'nwords'),
    ...        'print1': (print_and_return, 'format1'),
    ...        'print2': (print_and_return, 'format2'),
    ...        'print3': (print_and_return, 'format3')}

.. image:: images/optimize_dask1.svg
   :width: 65 %
   :alt: The original non-optimized Dask task graph.

Here we are counting the occurrence of the words ``'orange``, ``'apple'``, and
``'pear'`` in the list of words, formatting an output string reporting the
results, printing the output, and then returning the output string.

To perform the computation, we first remove unnecessary components from the
graph using the ``cull`` function and then pass the Dask graph and the desired
output keys to a scheduler ``get`` function:

.. code-block:: python

    >>> from dask.threaded import get
    >>> from dask.optimization import cull

    >>> outputs = ['print1', 'print2']
    >>> dsk1, dependencies = cull(dsk, outputs)  # remove unnecessary tasks from the graph

    >>> results = get(dsk1, outputs)
    word list has 2 occurrences of apple, out of 7 words
    word list has 2 occurrences of orange, out of 7 words

As can be seen above, the scheduler computed only the requested outputs
(``'print3'`` was never computed). This is because we called the
``dask.optimization.cull`` function, which removes the unnecessary tasks from
the graph.

Culling is part of the default optimization pass of almost all collections.
Often you want to call it somewhat early to reduce the amount of work done in
later steps:

.. code-block:: python

    >>> from dask.optimization import cull
    >>> outputs = ['print1', 'print2']
    >>> dsk1, dependencies = cull(dsk, outputs)

.. image:: images/optimize_dask2.svg
   :width: 50 %
   :alt: The Dask task graph after culling tasks for optimization.

Looking at the task graph above, there are multiple accesses to constants such
as ``'val1'`` or ``'val2'`` in the Dask graph. These can be inlined into the
tasks to improve efficiency using the ``inline`` function. For example:

.. code-block:: python

    >>> from dask.optimization import inline
    >>> dsk2 = inline(dsk1, dependencies=dependencies)
    >>> results = get(dsk2, outputs)
    word list has 2 occurrences of apple, out of 7 words
    word list has 2 occurrences of orange, out of 7 words

.. image:: images/optimize_dask3.svg
   :width: 40 %
   :alt: The Dask task graph after inlining for optimization.

Now we have two sets of *almost* linear task chains. The only link between them
is the word counting function. For cheap operations like this, the
serialization cost may be larger than the actual computation, so it may be
faster to do the computation more than once, rather than passing the results to
all nodes. To perform this function inlining, the ``inline_functions`` function
can be used:

.. code-block:: python

    >>> from dask.optimization import inline_functions
    >>> dsk3 = inline_functions(dsk2, outputs, [len, str.split],
    ...                         dependencies=dependencies)
    >>> results = get(dsk3, outputs)
    word list has 2 occurrences of apple, out of 7 words
    word list has 2 occurrences of orange, out of 7 words

.. image:: images/optimize_dask4.svg
   :width: 30 %
   :alt: The Dask task graph after inlining functions for optimization.

Now we have a set of purely linear tasks. We'd like to have the scheduler run
all of these on the same worker to reduce data serialization between workers.
One option is just to merge these linear chains into one big task using the
``fuse`` function:

.. code-block:: python

    >>> from dask.optimization import fuse
    >>> dsk4, dependencies = fuse(dsk3)
    >>> results = get(dsk4, outputs)
    word list has 2 occurrences of apple, out of 7 words
    word list has 2 occurrences of orange, out of 7 words

.. image:: images/optimize_dask5.svg
   :width: 30 %
   :alt: The Dask task graph after fusing tasks for optimization.


Putting it all together:

.. code-block:: python

    >>> def optimize_and_get(dsk, keys):
    ...     dsk1, deps = cull(dsk, keys)
    ...     dsk2 = inline(dsk1, dependencies=deps)
    ...     dsk3 = inline_functions(dsk2, keys, [len, str.split],
    ...                             dependencies=deps)
    ...     dsk4, deps = fuse(dsk3)
    ...     return get(dsk4, keys)

    >>> optimize_and_get(dsk, outputs)
    word list has 2 occurrences of apple, out of 7 words
    word list has 2 occurrences of orange, out of 7 words


In summary, the above operations accomplish the following:

1. Removed tasks unnecessary for the desired output using ``cull``
2. Inlined constants using ``inline``
3. Inlined cheap computations using ``inline_functions``, improving parallelism
4. Fused linear tasks together to ensure they run on the same worker using ``fuse``

As stated previously, these optimizations are already performed automatically
in the Dask collections. Users not working with custom graphs or computations
should rarely need to directly interact with them.

These are just a few of the optimizations provided in ``dask.optimization``. For
more information, see the API below.


Rewrite Rules
-------------

For context based optimizations, ``dask.rewrite`` provides functionality for
pattern matching and term rewriting. This is useful for replacing expensive
computations with equivalent, cheaper computations. For example, Dask Array
uses the rewrite functionality to replace series of array slicing operations
with a more efficient single slice.

The interface to the rewrite system consists of two classes:

1. ``RewriteRule(lhs, rhs, vars)``

    Given a left-hand-side (``lhs``), a right-hand-side (``rhs``), and a set of
    variables (``vars``), a rewrite rule declaratively encodes the following
    operation:

    ``lhs -> rhs if task matches lhs over variables``

2. ``RuleSet(*rules)``

    A collection of rewrite rules. The design of ``RuleSet`` class allows for
    efficient "many-to-one" pattern matching, meaning that there is minimal
    overhead for rewriting with multiple rules in a rule set.


Example
~~~~~~~

Here we create two rewrite rules expressing the following mathematical transformations:

1. ``a + a -> 2*a``
2. ``a * a -> a**2``

where ``'a'`` is a variable:

.. code-block:: python

    >>> from dask.rewrite import RewriteRule, RuleSet
    >>> from operator import add, mul, pow

    >>> variables = ('a',)

    >>> rule1 = RewriteRule((add, 'a', 'a'), (mul, 'a', 2), variables)

    >>> rule2 = RewriteRule((mul, 'a', 'a'), (pow, 'a', 2), variables)

    >>> rs = RuleSet(rule1, rule2)

The ``RewriteRule`` objects describe the desired transformations in a
declarative way, and the ``RuleSet`` builds an efficient automata for applying
that transformation. Rewriting can then be done using the ``rewrite`` method:

.. code-block:: python

    >>> rs.rewrite((add, 5, 5))
    (mul, 5, 2)

    >>> rs.rewrite((mul, 5, 5))
    (pow, 5, 2)

    >>> rs.rewrite((mul, (add, 3, 3), (add, 3, 3)))
    (pow, (mul, 3, 2), 2)

The whole task is traversed by default. If you only want to apply a transform
to the top-level of the task, you can pass in ``strategy='top_level'`` as shown:

.. code-block:: python

    # Transforms whole task
    >>> rs.rewrite((sum, [(add, 3, 3), (mul, 3, 3)]))
    (sum, [(mul, 3, 2), (pow, 3, 2)])

    # Only applies to top level, no transform occurs
    >>> rs.rewrite((sum, [(add, 3, 3), (mul, 3, 3)]), strategy='top_level')
    (sum, [(add, 3, 3), (mul, 3, 3)])

The rewriting system provides a powerful abstraction for transforming
computations at a task level. Again, for many users, directly interacting with
these transformations will be unnecessary.


Keyword Arguments
-----------------

Some optimizations take optional keyword arguments.  To pass keywords from the
compute call down to the right optimization, prepend the keyword with the name
of the optimization.  For example, to send a ``keys=`` keyword argument to the
``fuse`` optimization from a compute call, use the ``fuse_keys=`` keyword:

.. code-block:: python

   def fuse(dsk, keys=None):
       ...

   x.compute(fuse_keys=['x', 'y', 'z'])


Customizing Optimization
------------------------

Dask defines a default optimization strategy for each collection type (Array,
Bag, DataFrame, Delayed).  However, different applications may have different
needs.  To address this variability of needs, you can construct your own custom
optimization function and use it instead of the default.  An optimization
function takes in a task graph and list of desired keys and returns a new
task graph:

.. code-block:: python

   def my_optimize_function(dsk, keys):
       new_dsk = {...}
       return new_dsk

You can then register this optimization class against whichever collection type
you prefer and it will be used instead of the default scheme:

.. code-block:: python

   with dask.config.set(array_optimize=my_optimize_function):
       x, y = dask.compute(x, y)

You can register separate optimization functions for different collections, or
you can register ``None`` if you do not want particular types of collections to
be optimized:

.. code-block:: python

   with dask.config.set(array_optimize=my_optimize_function,
                        dataframe_optimize=None,
                        delayed_optimize=my_other_optimize_function):
       ...

You do not need to specify all collections.  Collections will default to their
standard optimization scheme (which is usually a good choice).


API
---

.. currentmodule:: dask.optimization

**Top level optimizations**

.. autosummary::
   cull
   fuse
   inline
   inline_functions

**Utility functions**

.. autosummary::
   functions_of

**Rewrite Rules**

.. currentmodule:: dask.rewrite

.. autosummary::
    RewriteRule
    RuleSet


Definitions
~~~~~~~~~~~

.. currentmodule:: dask.optimization

.. autofunction:: cull
.. autofunction:: fuse
.. autofunction:: inline
.. autofunction:: inline_functions

.. autofunction:: functions_of

.. currentmodule:: dask.rewrite

.. autofunction:: RewriteRule
.. autofunction:: RuleSet
.. _phases-of-computation:


Stages of Computation
=====================

This page describes all of the parts of computation, some common causes of
slowness, and how to effectively profile.  This is intended for more advanced
users who are encountering slowdowns on larger computations.

Graph Construction
------------------

Operations on Dask collections (array, dataframe, bag, delayed) build task
graphs.  These are dictionaries of Python functions that include an entry
every time some function needs to run on some chunk of data.  When these
dictionaries become large (millions of tasks) the overhead of constructing them
can become considerable.  Additionally the code that builds the graphs may
itself be inefficient.

Fortunately, this computation is all happening in normal Python right on your
own computer, so you can profile it just as you would any other Python code on
your computer using tools like the cProfile module, or the ``%prun`` or
``%snakeviz`` IPython magics.

Assuming that no obvious cause comes up when profiling, a common solution to
this problem is to reduce your graph size by increasing your chunk size if
possible, or manually batching many operations into fewer functions.

Graph Optimization
------------------

Just before you submit the graph to be executed, Dask sees if it can clean up
the graph a bit.  This helps to remove unnecessary work, and sometimes swaps
out more efficient operations.  As before though, if your graph is very large
(millions of tasks) then this can take some time.

Also as before, this is all happening in Python on your local machine.  You can
profile optimization separately from computation with the ``dask.optimize``
function.

.. code-block:: python

   # x, y = dask.compute(x, y)
   x, y = dask.optimize(x, y)

It's rare for people to change optimization.  It is rarely the main cause of
slowdown.


Graph serialization
-------------------

When you are using the distributed scheduler the graph must be sent to the
scheduler process, and from there to the workers.  To send this data it must
first be converted into bytes.  This serialization process can sometimes be
expensive either if the objects that you're passing around are very complex, or
if they are very large.

The easiest way to profile this is to profile the ``persist`` call with the
distributed scheduler.  This will include both the optimization phase above, as
well as the serialization and some of the communication phase below
(serialization is often the largest component).  Fortunately ``persist``
returns immediately after, not waiting for the computation to actually finish.

Most often the cause of long serialization times is placing large objects
like NumPy arrays or Pandas dataframes into your graph repeatedly.  Dask will
usually raise a warning when it notices this.  Often the best solution is to
read your data in as a task instead of include it directly, pre-scatter large
data, or wrap them in ``dask.delayed``.  Sometimes serialization is caused by
other issues with complex objects.  These tend to be very library specific, and
so it is hard to provide general guidelines for them.

Graph Communication
-------------------

The graph must then be communicated to the scheduler.  You can watch the
``/system`` tab of the dashboard to watch network communication to and from the
scheduler.  There is no good way to profile this.


Scheduling
----------

The scheduler now receives the graph, and must populate its internal data
structures to be able to efficiently schedule these tasks to the various
workers.

It is only after these data structures are populated that the dashboard will
show any activity.  All time between pressing ``compute/persist`` and seeing
activity is taken up in the stages above.

You can profile scheduling costs with the ``/profile-server`` page of the
dashboard.  However, this is rarely useful for users because, unless you're
willing to dive into the scheduling code, it is hard to act here.  Still, the
interested user may find the profile of the inter workings of the scheduler of
interest.

If scheduling is expensive then the best you can do is to reduce your graph
size, often by increasing chunk size.


Execution
---------

Finally your workers get sent some tasks and get to run them.  Your code runs
on a thread of a worker, and does whatever it was told to do.

Dask's Dashboard is a good tool to profile and investigate performance here,
particularly the ``/status`` and ``/profile`` pages.

Accelerating this phase is often up to the author of the tasks that you are
submitting.  This might be you if you are using custom code, or the NumPy or
Pandas developers.  We encourage you to consider efficient libraries like
Cython, Numba, or any other solution that is commonly used to accelerate Python
code.
Specification
=============

Dask is a specification to encode a graph -- specifically, a directed 
acyclic graph of tasks with data dependencies -- using ordinary Python data 
structures, namely dicts, tuples, functions, and arbitrary Python
values. 


Definitions
-----------

A **Dask graph** is a dictionary mapping **keys** to **computations**:

.. code-block:: python

   {'x': 1,
    'y': 2,
    'z': (add, 'x', 'y'),
    'w': (sum, ['x', 'y', 'z']),
    'v': [(sum, ['w', 'z']), 2]}

A **key** is any hashable value that is not a **task**:

.. code-block:: python

   'x'
   ('x', 2, 3)

A **task** is a tuple with a callable first element.  Tasks represent atomic
units of work meant to be run by a single worker.  Example: 

.. code-block:: python

   (add, 'x', 'y')

We represent a task as a tuple such that the *first element is a callable
function* (like ``add``), and the succeeding elements are *arguments* for that
function. An *argument* may be any valid **computation**.

A **computation** may be one of the following:

1.  Any **key** present in the Dask graph like ``'x'``
2.  Any other value like ``1``, to be interpreted literally
3.  A **task** like ``(inc, 'x')`` (see below)
4.  A list of **computations**, like ``[1, 'x', (inc, 'x')]``

So all of the following are valid **computations**:

.. code-block:: python

   np.array([...])
   (add, 1, 2)
   (add, 'x', 2)
   (add, (inc, 'x'), 2)
   (sum, [1, 2])
   (sum, ['x', (inc, 'x')])
   (np.dot, np.array([...]), np.array([...]))
   [(sum, ['x', 'y']), 'z']

To encode keyword arguments, we recommend the use of ``functools.partial`` or
``toolz.curry``.


What functions should expect
----------------------------

In cases like ``(add, 'x', 'y')``, functions like ``add`` receive concrete
values instead of keys.  A Dask scheduler replaces keys (like ``'x'`` and ``'y'``) with
their computed values (like ``1``, and ``2``) *before* calling the ``add`` function.


Entry Point - The ``get`` function
----------------------------------

The ``get`` function serves as entry point to computation for all
:doc:`schedulers <scheduler-overview>`.  This function gets the value
associated to the given key.  That key may refer to stored data, as is the case
with ``'x'``, or to a task, as is the case with ``'z'``.  In the latter case,
``get`` should perform all necessary computation to retrieve the computed
value.

.. _scheduler: scheduler-overview.rst

.. code-block:: python

   >>> from dask.threaded import get

   >>> from operator import add

   >>> dsk = {'x': 1,
   ...        'y': 2,
   ...        'z': (add, 'x', 'y'),
   ...        'w': (sum, ['x', 'y', 'z'])}

.. code-block:: python

   >>> get(dsk, 'x')
   1

   >>> get(dsk, 'z')
   3

   >>> get(dsk, 'w')
   6

Additionally, if given a ``list``, get should simultaneously acquire values for
multiple keys:

.. code-block:: python

   >>> get(dsk, ['x', 'y', 'z'])
   [1, 2, 3]

Because we accept lists of keys as keys, we support nested lists:

.. code-block:: python

   >>> get(dsk, [['x', 'y'], ['z', 'w']])
   [[1, 2], [3, 6]]

Internally ``get`` can be arbitrarily complex, calling out to distributed
computing, using caches, and so on.


Why use tuples
--------------

With ``(add, 'x', 'y')``, we wish to encode the result of calling ``add`` on
the values corresponding to the keys ``'x'`` and ``'y'``.

We intend the following meaning:

.. code-block:: python

   add('x', 'y')  # after x and y have been replaced

But this will err because Python executes the function immediately
before we know values for ``'x'`` and ``'y'``.

We delay the execution by moving the opening parenthesis one term to the left,
creating a tuple:

.. code::

    Before: add( 'x', 'y')
    After: (add, 'x', 'y')

This lets us store the desired computation as data that we can analyze using
other Python code, rather than cause immediate execution.

LISP users will identify this as an s-expression, or as a rudimentary form of
quoting.
Stats
=====

Dask Array implements a subset of the `scipy.stats`_ package.

Statistical Functions
---------------------

You can calculate various measures of an array including skewness, kurtosis, and arbitrary moments.

.. code-block:: python

   >>> from dask.array import stats
   >>> x = da.random.beta(1, 1, size=(1000,), chunks=10)
   >>> k, s, m = [stats.kurtosis(x), stats.skew(x), stats.moment(x, 5)]
   >>> dask.compute(k, s, m)
   (1.7612340817172787, -0.064073498030693302, -0.00054523780628304799)


Statistical Tests
-----------------

You can perform basic statistical tests on Dask arrays.
Each of these tests return a ``dask.delayed`` wrapping one of the scipy ``namedtuple``
results.


.. code-block:: python

   >>> a = da.random.uniform(size=(50,), chunks=(25,))
   >>> b = a + da.random.uniform(low=-0.15, high=0.15, size=(50,), chunks=(25,))
   >>> result = stats.ttest_rel(a, b)
   >>> result.compute()
   Ttest_relResult(statistic=-1.5102104380013242, pvalue=0.13741197274874514)

.. _scipy.stats: https://docs.scipy.org/doc/scipy-0.19.0/reference/stats.html
User Interfaces
===============

Dask supports several user interfaces:

-  High-Level
    -  :doc:`Arrays <array>`: Parallel NumPy
    -  :doc:`Bags <bag>`: Parallel lists
    -  :doc:`DataFrames <dataframe>`: Parallel Pandas
    -  `Machine Learning <https://ml.dask.org>`_ : Parallel Scikit-Learn
    -  Others from external projects, like `XArray <https://xarray.pydata.org>`_
-  Low-Level
    -  :doc:`Delayed <delayed>`: Parallel function evaluation
    -  :doc:`Futures <futures>`: Real-time parallel function evaluation

Each of these user interfaces employs the same underlying parallel computing
machinery, and so has the same scaling, diagnostics, resilience, and so on, but
each provides a different set of parallel algorithms and programming style.

This document helps you to decide which user interface best suits your needs,
and gives some general information that applies to all interfaces.
The pages linked above give more information about each interface in greater
depth.

High-Level Collections
----------------------

Many people who start using Dask are explicitly looking for a scalable version of
NumPy, Pandas, or Scikit-Learn.  For these situations, the starting point within
Dask is usually fairly clear.  If you want scalable NumPy arrays, then start with Dask
array; if you want scalable Pandas DataFrames, then start with Dask DataFrame, and so on.

These high-level interfaces copy the standard interface with slight variations.
These interfaces automatically parallelize over larger datasets for you for a
large subset of the API from the original project.

.. code-block:: python

   # Arrays
   import dask.array as da
   x = da.random.uniform(low=0, high=10, size=(10000, 10000),  # normal numpy code
                         chunks=(1000, 1000))  # break into chunks of size 1000x1000

   y = x + x.T - x.mean(axis=0)  # Use normal syntax for high level algorithms

   # DataFrames
   import dask.dataframe as dd
   df = dd.read_csv('2018-*-*.csv', parse_dates='timestamp',  # normal Pandas code
                    blocksize=64000000)  # break text into 64MB chunks

   s = df.groupby('name').balance.mean()  # Use normal syntax for high level algorithms

   # Bags / lists
   import dask.bag as db
   b = db.read_text('*.json').map(json.loads)
   total = (b.filter(lambda d: d['name'] == 'Alice')
             .map(lambda d: d['balance'])
             .sum())

It is important to remember that, while APIs may be similar, some differences do
exist.  Additionally, the performance of some algorithms may differ from their
in-memory counterparts due to the advantages and disadvantages of parallel
programming.  Some thought and attention is still required when using Dask.


Low-Level Interfaces
--------------------

Often when parallelizing existing code bases or building custom algorithms, you
run into code that is parallelizable, but isn't just a big DataFrame or array.
Consider the for-loopy code below:

.. code-block:: python

   results = []
   for a in A:
       for b in B:
           if a < b:
               c = f(a, b)
           else:
               c = g(a, b)
           results.append(c)

There is potential parallelism in this code (the many calls to ``f`` and ``g``
can be done in parallel), but it's not clear how to rewrite it into a big
array or DataFrame so that it can use a higher-level API.  Even if you could
rewrite it into one of these paradigms, it's not clear that this would be a
good idea.  Much of the meaning would likely be lost in translation, and this
process would become much more difficult for more complex systems.

Instead, Dask's lower-level APIs let you write parallel code one function call
at a time within the context of your existing for loops.  A common solution
here is to use :doc:`Dask delayed <delayed>` to wrap individual function calls
into a lazily constructed task graph:

.. code-block:: python

   import dask

   lazy_results = []
   for a in A:
       for b in B:
           if a < b:
               c = dask.delayed(f)(a, b)  # add lazy task
           else:
               c = dask.delayed(g)(a, b)  # add lazy task
           lazy_results.append(c)

   results = dask.compute(*lazy_results)  # compute all in parallel


Combining High- and Low-Level Interfaces
----------------------------------------

It is common to combine high- and low-level interfaces.
For example, you might use Dask array/bag/dataframe to load in data and do
initial pre-processing, then switch to Dask delayed for a custom algorithm that
is specific to your domain, then switch back to Dask array/dataframe to clean
up and store results.  Understanding both sets of user interfaces, and how
to switch between them, can be a productive combination.

.. code-block:: python

   # Convert to a list of delayed Pandas dataframes
   delayed_values = df.to_delayed()

   # Manipulate delayed values arbitrarily as you like

   # Convert many delayed Pandas DataFrames back to a single Dask DataFrame
   df = dd.from_delayed(delayed_values)


Laziness and Computing
----------------------

Most Dask user interfaces are *lazy*, meaning that they do not evaluate until
you explicitly ask for a result using the ``compute`` method:

.. code-block:: python

   # This array syntax doesn't cause computation
   y = x + x.T - x.mean(axis=0)

   # Trigger computation by explicitly calling the compute method
   y = y.compute()

If you have multiple results that you want to compute at the same time, use the
``dask.compute`` function.  This can share intermediate results and so be more
efficient:

.. code-block:: python

   # compute multiple results at the same time with the compute function
   min, max = dask.compute(y.min(), y.max())

Note that the ``compute()`` function returns in-memory results.  It converts
Dask DataFrames to Pandas DataFrames, Dask arrays to NumPy arrays, and Dask
bags to lists.  *You should only call compute on results that will fit
comfortably in memory*.  If your result does not fit in memory, then you might
consider writing it to disk instead.

.. code-block:: python

   # Write larger results out to disk rather than store them in memory
   my_dask_dataframe.to_parquet('myfile.parquet')
   my_dask_array.to_hdf5('myfile.hdf5')
   my_dask_bag.to_textfiles('myfile.*.txt')


Persist into Distributed Memory
-------------------------------

Alternatively, if you are on a cluster, then you may want to trigger a
computation and store the results in distributed memory.  In this case you do
not want to call ``compute``, which would create a single Pandas, NumPy, or
list result. Instead, you want to call ``persist``, which returns a new Dask
object that points to actively computing, or already computed results spread
around your cluster's memory.

.. code-block:: python

   # Compute returns an in-memory non-Dask object
   y = y.compute()

   # Persist returns an in-memory Dask object that uses distributed storage if available
   y = y.persist()

This is common to see after data loading an preprocessing steps, but before
rapid iteration, exploration, or complex algorithms.  For example, we might read
in a lot of data, filter down to a more manageable subset, and then persist
data into memory so that we can iterate quickly.

.. code-block:: python

   import dask.dataframe as dd
   df = dd.read_parquet('...')
   df = df[df.name == 'Alice']  # select important subset of data
   df = df.persist()  # trigger computation in the background

   # These are all relatively fast now that the relevant data is in memory
   df.groupby(df.id).balance.sum().compute()   # explore data quickly
   df.groupby(df.id).balance.mean().compute()  # explore data quickly
   df.id.nunique()                             # explore data quickly


Lazy vs Immediate
-----------------

As mentioned above, most Dask workloads are lazy, that is, they don't start any
work until you explicitly trigger them with a call to ``compute()``.
However, sometimes you *do* want to submit work as quickly as possible, track it
over time, submit new work or cancel work depending on partial results, and so
on.  This can be useful when tracking or responding to real-time events,
handling streaming data, or when building complex and adaptive algorithms.

For these situations, people typically turn to the :doc:`futures interface
<futures>` which is a low-level interface like Dask delayed, but operates
immediately rather than lazily.

Here is the same example with Dask delayed and Dask futures to illustrate the
difference.

Delayed: Lazy
~~~~~~~~~~~~~

.. code-block:: python

   @dask.delayed
   def inc(x):
       return x + 1

   @dask.delayed
   def add(x, y):
       return x + y

   a = inc(1)       # no work has happened yet
   b = inc(2)       # no work has happened yet
   c = add(a, b)    # no work has happened yet

   c = c.compute()  # This triggers all of the above computations


Futures: Immediate
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from dask.distributed import Client
   client = Client()

   def inc(x):
       return x + 1

   def add(x, y):
       return x + y

   a = client.submit(inc, 1)     # work starts immediately
   b = client.submit(inc, 2)     # work starts immediately
   c = client.submit(add, a, b)  # work starts immediately

   c = c.result()                # block until work finishes, then gather result

You can also trigger work with the high-level collections using the
``persist`` function.  This will cause work to happen in the background when
using the distributed scheduler.


Combining Interfaces
--------------------

There are established ways to combine the interfaces above:

1.  The high-level interfaces (array, bag, dataframe) have a ``to_delayed``
    method that can convert to a sequence (or grid) of Dask delayed objects

    .. code-block:: python

       delayeds = df.to_delayed()

2.  The high-level interfaces (array, bag, dataframe) have a ``from_delayed``
    method that can convert from either Delayed *or* Future objects

    .. code-block:: python

       df = dd.from_delayed(delayeds)
       df = dd.from_delayed(futures)

3.  The ``Client.compute`` method converts Delayed objects into Futures

    .. code-block:: python

       futures = client.compute(delayeds)

4.  The ``dask.distributed.futures_of`` function gathers futures from
    persisted collections

    .. code-block:: python

       from dask.distributed import futures_of

       df = df.persist()  # start computation in the background
       futures = futures_of(df)

5.  The Dask.delayed object converts Futures into delayed objects

    .. code-block:: python

       delayed_value = dask.delayed(future)

The approaches above should suffice to convert any interface into any other.
We often see some anti-patterns that do not work as well:

1.  Calling low-level APIs (delayed or futures) on high-level objects (like
    Dask arrays or DataFrames). This downgrades those objects to their NumPy or
    Pandas equivalents, which may not be desired.
    Often people are looking for APIs like ``dask.array.map_blocks`` or
    ``dask.dataframe.map_partitions`` instead.
2.  Calling ``compute()`` on Future objects.
    Often people want the ``.result()`` method instead.
3.  Calling NumPy/Pandas functions on high-level Dask objects or
    high-level Dask functions on NumPy/Pandas objects

Conclusion
----------

Most people who use Dask start with only one of the interfaces above but
eventually learn how to use a few interfaces together.  This helps them
leverage the sophisticated algorithms in the high-level interfaces while also
working around tricky problems with the low-level interfaces.

For more information, see the documentation for the particular user interfaces
below:

-  High Level
    -  :doc:`Arrays <array>`: Parallel NumPy
    -  :doc:`Bags <bag>`: Parallel lists
    -  :doc:`DataFrames <dataframe>`: Parallel Pandas
    -  `Machine Learning <https://ml.dask.org>`_ : Parallel Scikit-Learn
    -  Others from external projects, like `XArray <https://xarray.pydata.org>`_
-  Low Level
    -  :doc:`Delayed <delayed>`: Parallel function evaluation
    -  :doc:`Futures <futures>`: Real-time parallel function evaluation
SSH
===

It is easy to set up Dask on informally managed networks of machines using SSH.
This can be done manually using SSH and the
Dask :doc:`command line interface <deploying-cli>`,
or automatically using either the :class:`dask.distributed.SSHCluster` Python *cluster manager* or the
``dask-ssh`` command line tool. This document describes both of these options.

.. note::
   Before instaniating a ``SSHCluster`` it is recommended to configure keyless SSH
   for your local machine and other machines. For example, on a Mac to SSH into
   localhost (local machine) you need to ensure the Remote Login option is set in
   System Preferences -> Sharing. In addition, ``id_rsa.pub`` should be in
   ``authorized_keys`` for keyless login.

Python Interface
----------------

.. currentmodule:: dask.distributed

.. autofunction:: SSHCluster

Command Line
------------

The convenience script ``dask-ssh`` opens several SSH connections to your
target computers and initializes the network accordingly. You can
give it a list of hostnames or IP addresses::

   $ dask-ssh 192.168.0.1 192.168.0.2 192.168.0.3 192.168.0.4

Or you can use normal UNIX grouping::

   $ dask-ssh 192.168.0.{1,2,3,4}

Or you can specify a hostfile that includes a list of hosts::

   $ cat hostfile.txt
   192.168.0.1
   192.168.0.2
   192.168.0.3
   192.168.0.4

   $ dask-ssh --hostfile hostfile.txt

.. note::

   The command line documentation here may differ depending on your installed
   version. We recommend referring to the output of ``dask-ssh --help``.

.. click:: distributed.cli.dask_ssh:main
   :prog: dask-ssh
   :show-nested:
.. _dataframe.indexing:

Indexing into Dask DataFrames
=============================

Dask DataFrame supports some of Pandas' indexing behavior.

.. currentmodule:: dask.dataframe

.. autosummary::

   DataFrame.iloc
   DataFrame.loc

Label-based Indexing
--------------------

Just like Pandas, Dask DataFrame supports label-based indexing with the ``.loc``
accessor for selecting rows or columns, and ``__getitem__`` (square brackets)
for selecting just columns.

.. note::

   To select rows, the DataFrame's divisions must be known (see
   :ref:`dataframe.design` and :ref:`dataframe.performance` for more information.)

.. code-block:: python

   >>> import dask.dataframe as dd
   >>> import pandas as pd
   >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [3, 4, 5]},
   ...                   index=['a', 'b', 'c'])
   >>> ddf = dd.from_pandas(df, npartitions=2)
   >>> ddf
   Dask DataFrame Structure:
                      A      B
   npartitions=1
   a              int64  int64
   c                ...    ...
   Dask Name: from_pandas, 1 tasks

Selecting columns:

.. code-block:: python

   >>> ddf[['B', 'A']]
   Dask DataFrame Structure:
                      B      A
   npartitions=1
   a              int64  int64
   c                ...    ...
   Dask Name: getitem, 2 tasks


Selecting a single column reduces to a Dask Series:

.. code-block:: python

   >>> ddf['A']
   Dask Series Structure:
   npartitions=1
   a    int64
   c      ...
   Name: A, dtype: int64
   Dask Name: getitem, 2 tasks

Slicing rows and (optionally) columns with ``.loc``:

.. code-block:: python

   >>> ddf.loc[['b', 'c'], ['A']]
   Dask DataFrame Structure:
                      A
   npartitions=1
   b              int64
   c                ...
   Dask Name: loc, 2 tasks
   
   >>> ddf.loc[df["A"] > 1, ["B"]]
   Dask DataFrame Structure:
                      B
   npartitions=1
   a              int64
   c                ...
   Dask Name: try_loc, 2 tasks

   >>> ddf.loc[lambda df: df["A"] > 1, ["B"]]
   Dask DataFrame Structure:
                      B
   npartitions=1
   a              int64
   c                ...
   Dask Name: try_loc, 2 tasks

Dask DataFrame supports Pandas' `partial-string indexing <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#partial-string-indexing>`_:

.. code-block:: python

   >>> ts = dd.demo.make_timeseries()
   >>> ts
   Dask DataFrame Structure:
                      id    name        x        y
   npartitions=11
   2000-01-31      int64  object  float64  float64
   2000-02-29        ...     ...      ...      ...
   ...               ...     ...      ...      ...
   2000-11-30        ...     ...      ...      ...
   2000-12-31        ...     ...      ...      ...
   Dask Name: make-timeseries, 11 tasks

   >>> ts.loc['2000-02-12']
   Dask DataFrame Structure:
                                     id    name        x        y
   npartitions=1
   2000-02-12 00:00:00.000000000  int64  object  float64  float64
   2000-02-12 23:59:59.999999999    ...     ...      ...      ...
   Dask Name: loc, 12 tasks


Positional Indexing
-------------------

Dask DataFrame does not track the length of partitions, making positional
indexing with ``.iloc`` inefficient for selecting rows. :meth:`DataFrame.iloc`
only supports indexers where the row indexer is ``slice(None)`` (which ``:`` is
a shorthand for.)

.. code-block:: python

   >>> ddf.iloc[:, [1, 0]]
   Dask DataFrame Structure:
                      B      A
   npartitions=1
   a              int64  int64
   c                ...    ...
   Dask Name: iloc, 2 tasks

Trying to select specific rows with ``iloc`` will raise an exception:

.. code-block:: python

   >>> ddf.iloc[[0, 2], [1]]
   Traceback (most recent call last)
     File "<stdin>", line 1, in <module>
   ValueError: 'DataFrame.iloc' does not support slicing rows. The indexer must be a 2-tuple whose first item is 'slice(None)'.
Generalized Ufuncs
==================

`NumPy <https://www.numpy.org>`_ provides the concept of `generalized ufuncs <https://docs.scipy.org/doc/numpy/reference/c-api/generalized-ufuncs.html>`_. Generalized ufuncs are functions
that distinguish the various dimensions of passed arrays in the two classes loop dimensions
and core dimensions. To accomplish this, a `signature <https://docs.scipy.org/doc/numpy/reference/c-api/generalized-ufuncs.html#details-of-signature>`_ is specified for NumPy generalized ufuncs.

`Dask <https://dask.org/>`_ integrates interoperability with NumPy's generalized ufuncs
by adhering to respective `ufunc protocol <https://docs.scipy.org/doc/numpy/reference/arrays.classes.html#numpy.class.__array_ufunc__>`_, and provides a wrapper to make a Python function a generalized ufunc.


Usage
-----

NumPy Generalized UFuncs
~~~~~~~~~~~~~~~~~~~~~~~~
.. note::

    `NumPy <https://www.numpy.org>`_ generalized ufuncs are currently (v1.14.3 and below) stored in
    inside ``np.linalg._umath_linalg`` and might change in the future.


.. code-block:: python

    import dask.array as da
    import numpy as np

    x = da.random.normal(size=(3, 10, 10), chunks=(2, 10, 10))

    w, v = np.linalg._umath_linalg.eig(x, output_dtypes=(float, float))


Create Generalized UFuncs
~~~~~~~~~~~~~~~~~~~~~~~~~

It can be difficult to create your own GUFuncs without going into the CPython API.
However, the `Numba <https://numba.pydata.org>`_ project does provide a
nice implementation with their ``numba.guvectorize`` decorator.  See `Numba's
documentation
<https://numba.pydata.org/numba-doc/dev/user/vectorize.html#the-guvectorize-decorator>`_
for more information.

Wrap your own Python function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``gufunc`` can be used to make a Python function behave like a generalized ufunc:


.. code-block:: python

    x = da.random.normal(size=(10, 5), chunks=(2, 5))

    def foo(x):
        return np.mean(x, axis=-1)

    gufoo = da.gufunc(foo, signature="(i)->()", output_dtypes=float, vectorize=True)

    y = gufoo(x)


Instead of ``gufunc``, also the ``as_gufunc`` decorator can be used for convenience:


.. code-block:: python

    x = da.random.normal(size=(10, 5), chunks=(2, 5))

    @da.as_gufunc(signature="(i)->()", output_dtypes=float, vectorize=True)
    def gufoo(x):
        return np.mean(x, axis=-1)

    y = gufoo(x)


Disclaimer
----------
This experimental generalized ufunc integration is not complete:

* ``gufunc`` does not create a true generalized ufunc to be used with other input arrays besides Dask.
  I.e., at the moment, ``gufunc`` casts all input arguments to ``dask.array.Array``

* Inferring ``output_dtypes`` automatically is not implemented yet


API
---

.. currentmodule:: dask.array.gufunc

.. autosummary::
   apply_gufunc
   as_gufunc
   gufunc
Install Dask
============

You can install dask with ``conda``, with ``pip``, or by installing from source.

Conda
-----

Dask is installed by default in `Anaconda <https://www.anaconda.com/download/>`_.
You can update Dask using the `conda <https://www.anaconda.com/download/>`_ command::

   conda install dask

This installs Dask and **all** common dependencies, including Pandas and NumPy.
Dask packages are maintained both on the default channel and on `conda-forge <https://conda-forge.github.io/>`_.
Optionally, you can obtain a minimal Dask installation using the following command::

   conda install dask-core

This will install a minimal set of dependencies required to run Dask similar to (but not exactly the same as) ``python -m pip install dask`` below.

Pip
---

You can install everything required for most common uses of Dask (arrays,
dataframes, ...)  This installs both Dask and dependencies like NumPy, Pandas,
and so on that are necessary for different workloads.  This is often the right
choice for Dask users::

   python -m pip install "dask[complete]"    # Install everything

You can also install only the Dask library.  Modules like ``dask.array``,
``dask.dataframe``, or ``dask.distributed`` won't work until you also install NumPy,
Pandas, or Tornado, respectively.  This is common for downstream library
maintainers::

   python -m pip install dask                # Install only core parts of dask

We also maintain other dependency sets for different subsets of functionality::

   python -m pip install "dask[array]"       # Install requirements for dask array
   python -m pip install "dask[dataframe]"   # Install requirements for dask dataframe
   python -m pip install "dask[diagnostics]" # Install requirements for dask diagnostics
   python -m pip install "dask[distributed]" # Install requirements for distributed dask

We have these options so that users of the lightweight core Dask scheduler
aren't required to download the more exotic dependencies of the collections
(Numpy, Pandas, Tornado, etc.).


Install from Source
-------------------

To install Dask from source, clone the repository from `github
<https://github.com/dask/dask>`_::

    git clone https://github.com/dask/dask.git
    cd dask
    python -m pip install .

You can also install all dependencies as well::

    python -m pip install ".[complete]"

You can view the list of all dependencies within the ``extras_require`` field
of ``setup.py``.


Or do a developer install by using the ``-e`` flag::

    python -m pip install -e .

Anaconda
--------

Dask is included by default in the `Anaconda distribution <https://www.anaconda.com/download>`_.

Optional dependencies
---------------------

Specific functionality in Dask may require additional optional dependencies.
For example, reading from Amazon S3 requires ``s3fs``.
These optional dependencies and their minimum supported versions are listed below.

+---------------+----------+--------------------------------------------------------------+
| Dependency    | Version  |                          Description                         |
+===============+==========+==============================================================+
|     bokeh     | >=2.1.1  |                Visualizing dask diagnostics                  |
+---------------+----------+--------------------------------------------------------------+
|   cityhash    |          |                  Faster hashing of arrays                    |
+---------------+----------+--------------------------------------------------------------+
|  distributed  | >=2.0    |               Distributed computing in Python                |
+---------------+----------+--------------------------------------------------------------+
|  fastparquet  |          |         Storing and reading data from parquet files          |
+---------------+----------+--------------------------------------------------------------+
|     gcsfs     | >=0.4.0  |        File-system interface to Google Cloud Storage         |
+---------------+----------+--------------------------------------------------------------+
|   murmurhash  |          |                   Faster hashing of arrays                   |
+---------------+----------+--------------------------------------------------------------+
|     numpy     | >=1.18   |                   Required for dask.array                    |
+---------------+----------+--------------------------------------------------------------+
|     pandas    | >=1.0    |                  Required for dask.dataframe                 |
+---------------+----------+--------------------------------------------------------------+
|     psutil    |          |             Enables a more accurate CPU count                |
+---------------+----------+--------------------------------------------------------------+
|     pyarrow   | >=1.0    |               Python library for Apache Arrow                |
+---------------+----------+--------------------------------------------------------------+
|     s3fs      | >=0.4.0  |                    Reading from Amazon S3                    |
+---------------+----------+--------------------------------------------------------------+
|     scipy     |          |                  Required for dask.array.stats               |
+---------------+----------+--------------------------------------------------------------+
|   sqlalchemy  |          |            Writing and reading from SQL databases            |
+---------------+----------+--------------------------------------------------------------+
|    cytoolz*   | >=0.8.2  | Utility functions for iterators, functions, and dictionaries |
+---------------+----------+--------------------------------------------------------------+
|    xxhash     |          |                  Faster hashing of arrays                    |
+---------------+----------+--------------------------------------------------------------+

\* Note that ``toolz`` is a mandatory dependency but it can be transparently replaced with
``cytoolz``.


Test
----

Test Dask with ``py.test``::

    cd dask
    py.test dask

Please be aware that installing Dask naively may not install all
requirements by default. Please read the ``pip`` section above which discusses
requirements.  You may choose to install the ``dask[complete]`` version which includes
all dependencies for all collections.  Alternatively, you may choose to test
only certain submodules depending on the libraries within your environment.
For example, to test only Dask core and Dask array we would run tests as
follows::

    py.test dask/tests dask/array/tests
Create and Store Dask DataFrames
================================

Dask can create DataFrames from various data storage formats like CSV, HDF,
Apache Parquet, and others.  For most formats, this data can live on various
storage systems including local disk, network file systems (NFS), the Hadoop
File System (HDFS), and Amazon's S3 (excepting HDF, which is only available on
POSIX like file systems).

See the :doc:`DataFrame overview page <dataframe>` for an in depth
discussion of ``dask.dataframe`` scope, use, and limitations.

API
---

The following functions provide access to convert between Dask DataFrames,
file formats, and other Dask or Python collections.

.. currentmodule:: dask.dataframe

File Formats:

.. autosummary::
    read_csv
    read_parquet
    read_hdf
    read_orc
    read_json
    read_sql_table
    read_sql_query
    read_sql
    read_table
    read_fwf
    from_bcolz
    from_array
    to_csv
    to_parquet
    to_hdf
    to_sql

Dask Collections:

.. autosummary::
    from_delayed
    from_dask_array
    dask.bag.core.Bag.to_dataframe
    DataFrame.to_delayed
    to_records
    to_bag

Pandas:

.. autosummary::
    from_pandas

Creating
--------

Reading from various locations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For text, CSV, and Apache Parquet formats, data can come from local disk,
the Hadoop File System, S3FS, or other sources, by prepending the filenames with
a protocol:

.. code-block:: python

   >>> df = dd.read_csv('my-data-*.csv')
   >>> df = dd.read_csv('hdfs:///path/to/my-data-*.csv')
   >>> df = dd.read_csv('s3://bucket-name/my-data-*.csv')

For remote systems like HDFS, S3 or GS credentials may be an issue.  Usually, these
are handled by configuration files on disk (such as a ``.boto`` file for S3),
but in some cases you may want to pass storage-specific options through to the
storage backend.  You can do this with the ``storage_options=`` keyword:

.. code-block:: python

   >>> df = dd.read_csv('s3://bucket-name/my-data-*.csv',
   ...                  storage_options={'anon': True})
   >>> df = dd.read_parquet('gs://dask-nyc-taxi/yellowtrip.parquet',
   ...                      storage_options={'token': 'anon'})

Dask Delayed
~~~~~~~~~~~~

For more complex situations not covered by the functions above, you may want to
use :doc:`dask.delayed<delayed>`, which lets you construct Dask DataFrames out
of arbitrary Python function calls that load DataFrames.  This can allow you to
handle new formats easily or bake in particular logic around loading data if,
for example, your data is stored with some special format.

See :doc:`documentation on using dask.delayed with
collections<delayed-collections>` or an `example notebook
<https://gist.github.com/mrocklin/e7b7b3a65f2835cda813096332ec73ca>`_ showing
how to create a Dask DataFrame from a nested directory structure of Feather
files (as a stand in for any custom file format).

Dask delayed is particularly useful when simple ``map`` operations aren't
sufficient to capture the complexity of your data layout.

From Raw Dask Graphs
~~~~~~~~~~~~~~~~~~~~

This section is mainly for developers wishing to extend ``dask.dataframe``.  It
discusses internal API not normally needed by users.  Everything below can be
done just as effectively with :doc:`dask.delayed<delayed>`  described
just above.  You should never need to create a DataFrame object by hand.

To construct a DataFrame manually from a dask graph you need the following
information:

1.  Dask: a Dask graph with keys like ``{(name, 0): ..., (name, 1): ...}`` as
    well as any other tasks on which those tasks depend.  The tasks
    corresponding to ``(name, i)`` should produce ``pandas.DataFrame`` objects
    that correspond to the columns and divisions information discussed below
2.  Name: the special name used above
3.  Meta: an empty pandas DataFrame with names, dtypes and index matching
    the expected output. Can also be a list of tuples where each tuple defines
    a ``(name, dtype)`` pair referring to one column.
4.  Divisions: a list of index values that separate the different partitions.
    Alternatively, if you don't know the divisions (this is common), you can
    provide a list of ``[None, None, None, ...]`` with as many partitions as
    you have plus one.  For more information, see the Partitions section in the
    :doc:`DataFrame documentation <dataframe>`

As an example, we build a DataFrame manually that reads several CSV files that
have a datetime index separated by day.  Note that you should **never** do this.
The ``dd.read_csv`` function does this for you:

.. code-block:: Python

   dsk = {('mydf', 0): (pd.read_csv, 'data/2000-01-01.csv'),
          ('mydf', 1): (pd.read_csv, 'data/2000-01-02.csv'),
          ('mydf', 2): (pd.read_csv, 'data/2000-01-03.csv')}
   name = 'mydf'
   meta = [('price', float), ('name', str), ('id', int)]
   divisions = [Timestamp('2000-01-01 00:00:00'),
                Timestamp('2000-01-02 00:00:00'),
                Timestamp('2000-01-03 00:00:00'),
                Timestamp('2000-01-03 23:59:59')]

   df = dd.DataFrame(dsk, name, meta, divisions)

Storing
-------

Writing to remote locations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask can write to a variety of data stores including cloud object stores.
For example, you can write a ``dask.dataframe`` to an Azure storage blob as:

.. code-block:: python

   >>> d = {'col1': [1, 2, 3, 4], 'col2': [5, 6, 7, 8]}
   >>> df = dd.from_pandas(pd.DataFrame(data=d), npartitions=2)
   >>> dd.to_parquet(df=df,
   ...               path='abfs://CONTAINER/FILE.parquet'
   ...               storage_options={'account_name': 'ACCOUNT_NAME',
   ...                                'account_key': 'ACCOUNT_KEY'}

See the :doc:`how to connect to remote data <how-to/connect-to-remote-data>`
for more information.
Talks & Tutorials
=================

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/nnndxbr_Xq4"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

Dask Tutorial
-------------
`Dask Tutorial <https://tutorial.dask.org>`__ provides an overview of Dask and is typically delivered in 3 hours.
See `Parallel and Distributed Computing in Python with Dask <https://www.youtube.com/watch?v=EybGGLbLipI>`__ for the
latest Dask Tutorial recording from SciPy 2020.

Dask Slides
-----------
`Dask Slides <https://dask.org/slides>`__ provide a quick overview of the motivation for Dask.

Dask YouTube channel
--------------------
You can find lots of videos about Dask on the `Dask YouTube channel <https://www.youtube.com/c/dask-dev>`__

.. contents:: :local:

Presentations
-------------

* Dask Summit 2021

  * `Keynotes <https://www.youtube.com/playlist?list=PLJ0vO2F_f6OBymP5LtgOC6W4pxd9Mw3cE>`__ 
  * `Workshops and Tutorials <https://www.youtube.com/playlist?list=PLJ0vO2F_f6OBD1_iNeT1f7cpRoYwAuMPy>`__
  * `Talks <https://www.youtube.com/playlist?list=PLJ0vO2F_f6OBcisTDubrdEsQAhigkayjE>`__

* PyCon US 2021

  * `Tutorial: Hacking Dask: Diving into Dask's Internals <https://www.youtube.com/watch?v=LQrgDhN-XOo>`__  (`materials <https://github.com/jrbourbeau/hacking-dask>`__)
  * `Dask-SQL: Empowering Pythonistas for Scalable End-to-End Data Engineering <https://www.youtube.com/watch?v=z7xKikaScxg>`__


* BlazingSQL Webinars, May 2021

  * `Intro to distributed computing on GPUs with Dask in Python <https://www.youtube.com/watch?v=py1YPs6s6so>`__ (`materials <https://gist.github.com/jacobtomlinson/6f16abb716f50f81a6687bd67efd2f61>`__)

* PyData DC, August 2021

  * `Inside Dask <https://www.youtube.com/watch?v=X95WO41abXo>`__ (`materials <https://github.com/jsignell/inside-dask>`__)

* PyCon US 2020

  * `Deploying Python at Scale with Dask <https://www.youtube.com/watch?v=deX0GlW4uew>`__

* PyCon Australia 2020

  * `dask-image: distributed image processing for large data <https://www.youtube.com/watch?v=MpjgzNeISeI>`__

* PyCon Korea 2019, August 2019

  * `Adapting from Spark to Dask: what to expect (18 minutes)
    <https://www.youtube.com/watch?v=tx7qTHSlHKw>`__

* SciPy 2019, July 2019

  * `Refactoring the SciPy Ecosystem for Heterogeneous Computing (29 minutes)
    <https://www.youtube.com/watch?v=Q0DsdiY-jiw>`__
  * `Renewable Power Forecast Generation with Dask & Visualization with Bokeh (31 minutes)
    <https://www.youtube.com/watch?v=tYGcicSruck>`__
  * `Efficient Atmospheric Analogue Selection with Xarray and Dask (18 minutes)
    <https://www.youtube.com/watch?v=gdHiGsGUh3o>`__
  * `Better and Faster Hyper Parameter Optimization with Dask (27 minutes)
    <https://www.youtube.com/watch?v=x67K9FiPFBQ>`__
  * `Dask image:A Library for Distributed Image Processing (22 minutes)
    <https://www.youtube.com/watch?v=XGUS174vvLs>`__

* EuroPython 2019, July 2019

  * `Distributed Multi-GPU Computing with Dask, CuPy and RAPIDS (29 minutes)
    <https://www.youtube.com/watch?v=en2zdTT-Vwk>`__

* SciPy 2018, July 2018

  * `Scalable Machine Learning with Dask (30 minutes)
    <https://www.youtube.com/watch?v=ccfsbuqsjgI>`__

* PyCon 2018, May 2018

  * `Democratizing Distributed Computing with Dask and JupyterHub (32 minutes)
    <https://www.youtube.com/watch?v=Iq72dt1gO9c>`__

* AMS & ESIP, January 2018

  * `Pangeo quick demo: Dask, XArray, Zarr on the cloud with JupyterHub (3 minutes)
    <https://www.youtube.com/watch?v=rSOJKbfNBNk>`__
  * `Pangeo talk: An open-source big data science platform with Dask, XArray, Zarr on the cloud with JupyterHub (43 minutes)
    <https://www.youtube.com/watch?v=mDrjGxaXQT4>`__

* PYCON.DE 2017, November 2017

  * `Dask: Parallelism in Python (1 hour, 2 minutes)
    <https://www.youtube.com/watch?v=rZlshXJydgQ>`__

* PYCON 2017, May 2017

  * `Dask: A Pythonic Distributed Data Science Framework (46 minutes)
    <https://www.youtube.com/watch?v=RA_2qdipVng>`__

* PLOTCON 2016, December 2016

  * `Visualizing Distributed Computations with Dask and Bokeh (33 minutes)
    <https://www.youtube.com/watch?v=FTJwDeXkggU>`__

* PyData DC, October 2016

  * `Using Dask for Parallel Computing in Python (44 minutes)
    <https://www.youtube.com/watch?v=s4ChP7tc3tA>`__

* SciPy 2016, July 2016

  * `Dask Parallel and Distributed Computing (28 minutes)
    <https://www.youtube.com/watch?v=PAGjm4BMKlk>`__

* PyData NYC, December 2015

  * `Dask Parallelizing NumPy and Pandas through Task Scheduling (33 minutes)
    <https://www.youtube.com/watch?v=mHd8AI8GQhQ>`__

* PyData Seattle, August 2015

  * `Dask: out of core arrays with task scheduling (1 hour, 50 minutes)
    <https://www.youtube.com/watch?v=ieW3G7ZzRZ0>`__

* SciPy 2015, July 2015

  * `Dask Out of core NumPy:Pandas through Task Scheduling (16 minutes)
    <https://www.youtube.com/watch?v=1kkFZ4P-XHg>`__
:orphan:

Dask Cheat Sheet
================

The 300KB pdf :download:`Dask cheat sheet <daskcheatsheet.pdf>`
is a single page summary about using Dask.
It is commonly distributed at conferences and trade shows.
DataFrame
=========

.. toctree::
   :maxdepth: 1
   :hidden:

   dataframe-create.rst
   dataframe-best-practices.rst
   dataframe-design.rst
   dataframe-groupby.rst
   dataframe-joins.rst
   dataframe-indexing.rst
   dataframe-categoricals.rst
   dataframe-extend.rst
   dataframe-sql.rst

A Dask DataFrame is a large parallel DataFrame composed of many smaller Pandas
DataFrames, split along the index.  These Pandas DataFrames may live on disk
for larger-than-memory computing on a single machine, or on many different
machines in a cluster.  One Dask DataFrame operation triggers many operations
on the constituent Pandas DataFrames.

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/AT2XtFehFSQ"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>

Examples
--------

Visit https://examples.dask.org/dataframe.html to see and run examples using
Dask DataFrame.

Design
------

.. image:: images/dask-dataframe.svg
   :alt: Column of four squares collectively labeled as a Dask DataFrame with a single constituent square labeled as a Pandas DataFrame.
   :width: 35%
   :align: right

Dask DataFrames coordinate many Pandas DataFrames/Series arranged along the
index.  A Dask DataFrame is partitioned *row-wise*, grouping rows by index value
for efficiency.  These Pandas objects may live on disk or on other machines.

|

Dask DataFrame copies the Pandas API
------------------------------------

Because the ``dask.dataframe`` application programming interface (API) is a
subset of the Pandas API, it should be familiar to Pandas users.  There are some
slight alterations due to the parallel nature of Dask:

.. code-block:: python

   >>> import dask.dataframe as dd
   >>> df = dd.read_csv('2014-*.csv')
   >>> df.head()
      x  y
   0  1  a
   1  2  b
   2  3  c
   3  4  a
   4  5  b
   5  6  c

   >>> df2 = df[df.y == 'a'].x + 1

As with all Dask collections, one triggers computation by calling the
``.compute()`` method:

.. code-block:: python

   >>> df2.compute()
   0    2
   3    5
   Name: x, dtype: int64


Common Uses and Anti-Uses
-------------------------

Dask DataFrame is used in situations where Pandas is commonly needed, usually when
Pandas fails due to data size or computation speed:

-  Manipulating large datasets, even when those datasets don't fit in memory
-  Accelerating long computations by using many cores
-  Distributed computing on large datasets with standard Pandas operations like
   groupby, join, and time series computations

Dask DataFrame may not be the best choice in the following situations:

*  If your dataset fits comfortably into RAM on your laptop, then you may be
   better off just using Pandas.  There may be simpler ways to improve
   performance than through parallelism
*  If your dataset doesn't fit neatly into the Pandas tabular model, then you
   might find more use in :doc:`dask.bag <bag>` or :doc:`dask.array <array>`
*  If you need functions that are not implemented in Dask DataFrame, then you
   might want to look at :doc:`dask.delayed <delayed>` which offers more
   flexibility
*  If you need a proper database with all that databases offer you might prefer
   something like Postgres_

.. _Pandas: https://pandas.pydata.org/
.. _Postgres: https://www.postgresql.org/


Scope
-----

Dask DataFrame covers a well-used portion of the Pandas API.
The following class of computations works well:

* Trivially parallelizable operations (fast):
    *  Element-wise operations:  ``df.x + df.y``, ``df * df``
    *  Row-wise selections:  ``df[df.x > 0]``
    *  Loc:  ``df.loc[4.0:10.5]``
    *  Common aggregations:  ``df.x.max()``, ``df.max()``
    *  Is in:  ``df[df.x.isin([1, 2, 3])]``
    *  Date time/string accessors:  ``df.timestamp.month``
* Cleverly parallelizable operations (fast):
    *  groupby-aggregate (with common aggregations): ``df.groupby(df.x).y.max()``,
       ``df.groupby('x').max()``
    *  groupby-apply on index: ``df.groupby(['idx', 'x']).apply(myfunc)``, where
       ``idx`` is the index level name
    *  value_counts:  ``df.x.value_counts()``
    *  Drop duplicates:  ``df.x.drop_duplicates()``
    *  Join on index:  ``dd.merge(df1, df2, left_index=True, right_index=True)``
       or ``dd.merge(df1, df2, on=['idx', 'x'])`` where ``idx`` is the index
       name for both ``df1`` and ``df2``
    *  Join with Pandas DataFrames: ``dd.merge(df1, df2, on='id')``
    *  Element-wise operations with different partitions / divisions: ``df1.x + df2.y``
    *  Date time resampling: ``df.resample(...)``
    *  Rolling averages:  ``df.rolling(...)``
    *  Pearson's correlation: ``df[['col1', 'col2']].corr()``
* Operations requiring a shuffle (slow-ish, unless on index)
    *  Set index:  ``df.set_index(df.x)``
    *  groupby-apply not on index (with anything):  ``df.groupby(df.x).apply(myfunc)``
    *  Join not on the index:  ``dd.merge(df1, df2, on='name')``

However, Dask DataFrame does not implement the entire Pandas interface.  Users
expecting this will be disappointed.  Notably, Dask DataFrame has the following
limitations:

1.  Setting a new index from an unsorted column is expensive
2.  Many operations like groupby-apply and join on unsorted columns require
    setting the index, which as mentioned above, is expensive
3.  The Pandas API is very large.  Dask DataFrame does not attempt to implement
    many Pandas features or any of the more exotic data structures like NDFrames
4.  Operations that were slow on Pandas, like iterating through row-by-row,
    remain slow on Dask DataFrame

See :doc:`DataFrame API documentation<dataframe-api>` for a more extensive list.


Execution
---------

By default, Dask DataFrame uses the multi-threaded scheduler.
This exposes some parallelism when Pandas or the underlying NumPy operations
release the global interpreter lock (GIL).  Generally, Pandas is more GIL
bound than NumPy, so multi-core speed-ups are not as pronounced for
Dask DataFrame as they are for Dask Array.  This is changing, and
the Pandas development team is actively working on releasing the GIL.

When dealing with text data, you may see speedups by switching to the
:doc:`distributed scheduler <how-to/deploy-dask/single-distributed>` either on a cluster or
single machine.
LinearOperator
==============

Dask Array implements the SciPy LinearOperator_ interface and it can be used
with any SciPy algorithm depending on that interface.

Example
-------

.. code-block:: python

   import dask.array as da
   x = da.random.random(size=(10000, 10000), chunks=(1000, 1000))

   from scipy.sparse.linalg.interface import MatrixLinearOperator
   A = MatrixLinearOperator(x)

   import numpy as np
   b = np.random.random(10000)

   from scipy.sparse.linalg import gmres
   x = gmres(A, b)

*Disclaimer: This is just a toy example and not necessarily the best way to
solve this problem for this data.*


.. _LinearOperator: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.LinearOperator.html
Create Dask Bags
================

There are several ways to create Dask bags around your data:

``db.from_sequence``
--------------------

You can create a bag from an existing Python iterable:

.. code-block:: python

   >>> import dask.bag as db
   >>> b = db.from_sequence([1, 2, 3, 4, 5, 6])

You can control the number of partitions into which this data is binned:

.. code-block:: python

   >>> b = db.from_sequence([1, 2, 3, 4, 5, 6], npartitions=2)

This controls the granularity of the parallelism that you expose.  By default,
Dask will try to partition your data into about 100 partitions.

IMPORTANT: do not load your data into Python and then load that data into a
Dask bag.  Instead, use Dask Bag to load your data.  This
parallelizes the loading step and reduces inter-worker communication:

.. code-block:: python

   >>> b = db.from_sequence(['1.dat', '2.dat', ...]).map(load_from_filename)


``db.read_text``
----------------

Dask Bag can load data directly from text files.  You can pass either a 
single file name, a list of file names, or a globstring.  The resulting 
bag will have one item per line and one file per partition:

.. code-block:: python

   >>> b = db.read_text('myfile.txt')
   >>> b = db.read_text(['myfile.1.txt', 'myfile.2.txt', ...])
   >>> b = db.read_text('myfile.*.txt')

This handles standard compression libraries like ``gzip``, ``bz2``, ``xz``, or
any easily installed compression library that has a file-like object.
Compression will be inferred by the file name extension, or by using the
``compression='gzip'`` keyword:

.. code-block:: python

   >>> b = db.read_text('myfile.*.txt.gz')

The resulting items in the bag are strings. If you have encoded data like
line-delimited JSON, then you may want to map a decoding or load function across
the bag:

.. code-block:: python

   >>> import json
   >>> b = db.read_text('myfile.*.json').map(json.loads)

Or do string munging tasks.  For convenience, there is a string namespace
attached directly to bags with ``.str.methodname``:

.. code-block:: python

   >>> b = db.read_text('myfile.*.csv').str.strip().str.split(',')


``db.read_avro``
----------------

Dask Bag can read binary files in the `Avro`_ format if `fastavro`_ is installed.
A bag can be made from one or more files, with optional chunking within files. 
The resulting bag will have one item per Avro record, which will be a dictionary 
of the form given by the Avro schema.  There will be at least one partition per 
input file:

.. code-block:: python

   >>> b = db.read_avro('datafile.avro')
   >>> b = db.read_avro('data.*.avro')

.. _Avro: https://avro.apache.org/docs/1.8.2/
.. _fastavro: https://fastavro.readthedocs.io

By default, Dask will split data files into chunks of approximately ``blocksize`` 
bytes in size.  The actual blocks you would get depend on the internal blocking 
of the file.

For files that are compressed after creation (this is not the same as the internal 
"codec" used by Avro), no chunking should be used, and there will be exactly one 
partition per file:

.. code-block:: python

   > b = bd.read_avro('compressed.*.avro.gz', blocksize=None, compression='gzip')


``db.from_delayed``
-------------------

You can construct a Dask bag from :doc:`dask.delayed <delayed>` values using the 
``db.from_delayed`` function.  For more information, see
:doc:`documentation on using dask.delayed with collections <delayed-collections>`.


Store Dask Bags
===============

In Memory
---------

You can convert a Dask bag to a list or Python iterable by calling ``compute()`` 
or by converting the object into a list:

.. code-block:: python

   >>> result = b.compute()
   or
   >>> result = list(b)

To Text Files
-------------

You can convert a Dask bag into a sequence of files on disk by calling the
``.to_textfiles()`` method:

.. autofunction:: dask.bag.core.to_textfiles

To Avro
-------

Dask bags can be written directly to Avro binary format using `fastavro`_. One file
will be written per bag partition. This requires the user to provide a fully-specified
schema dictionary (see the docstring of the ``.to_avro()`` method).

.. autofunction:: dask.bag.avro.to_avro

To DataFrames
-------------

You can convert a Dask bag into a :doc:`Dask DataFrame<dataframe>` and use
those storage solutions.

.. automethod:: dask.bag.core.Bag.to_dataframe


To Delayed Values
-----------------

You can convert a Dask bag into a list of :doc:`Dask delayed values<delayed>`
and custom storage solutions from there.

.. automethod:: dask.bag.core.Bag.to_delayed
Changelog
=========

.. _v2022.01.1:

2022.01.1
---------

Released on January 28, 2022

New Features
^^^^^^^^^^^^
- Add ``dask.dataframe.series.view()`` (:pr:`8533`) `Pavithra Eswaramoorthy`_

Enhancements
^^^^^^^^^^^^
- Update ``tz`` for ``fastparquet`` + ``pandas`` 1.4.0 (:pr:`8626`) `Martin Durant`_
- Cleaning up misc tests for ``pandas`` compat (:pr:`8623`) `Julia Signell`_
- Moving to ``SQLAlchemy >= 1.4`` (:pr:`8158`) `McToel`_
- Pandas compat: Filter sparse warnings (:pr:`8621`) `Julia Signell`_
- Fail if ``meta`` is not a ``pandas`` object (:pr:`8563`) `Julia Signell`_
- Use ``fsspec.parquet`` module for better remote-storage ``read_parquet`` performance (:pr:`8339`) `Richard (Rick) Zamora`_
- Move DataFrame ACA aggregations to HLG (:pr:`8468`) `Richard (Rick) Zamora`_
- Add optional information about originating function call in ``DataFrameIOLayer`` (:pr:`8453`) `Richard (Rick) Zamora`_
- Blockwise array creation redux (:pr:`7417`) `Ian Rose`_
- Refactor config default search path retrieval (:pr:`8573`) `James Bourbeau`_
- Add ``optimize_graph`` flag to ``Bag.to_dataframe`` function (:pr:`8486`) `Maxim Lippeveld`_
- Make sure that delayed output operations still return lists of paths (:pr:`8498`) `Julia Signell`_
- Pandas compat: Fix ``to_frame`` ``name`` to not pass ``None`` (:pr:`8554`) `Julia Signell`_
- Pandas compat: Fix ``axis=None`` warning (:pr:`8555`) `Julia Signell`_
- Expand Dask YAML config search directories (:pr:`8531`) `abergou`_

Bug Fixes
^^^^^^^^^
- Fix ``groupby.cumsum`` with series grouped by index (:pr:`8588`) `Julia Signell`_
- Fix ``derived_from`` for ``pandas`` methods (:pr:`8612`) `Thomas J. Fan`_
- Enforce boolean ``ascending`` for ``sort_values`` (:pr:`8440`) `Charles Blackmon-Luca`_
- Fix parsing of ``__setitem__`` indices (:pr:`8601`) `David Hassell`_
- Avoid divide by zero in slicing  (:pr:`8597`) `Doug Davis`_

Deprecations
^^^^^^^^^^^^
- Downgrade ``meta`` error in (:pr:`8563`) to warning (:pr:`8628`) `Julia Signell`_
- Pandas compat: Deprecate ``append`` when ``pandas >= 1.4.0`` (:pr:`8617`) `Julia Signell`_

Documentation
^^^^^^^^^^^^^
- Replace outdated ``columns`` argument with ``meta`` in DataFrame constructor (:pr:`8614`) `kori73`_
- Refactor deploying docs (:pr:`8602`) `Jacob Tomlinson`_

Maintenance
^^^^^^^^^^^
- Pin ``coverage`` in CI (:pr:`8631`) `James Bourbeau`_
- Move ``cached_cumsum`` imports to be from ``dask.utils`` (:pr:`8606`) `James Bourbeau`_
- Update gpuCI ``RAPIDS_VER`` to ``22.04`` (:pr:`8600`)
- Update cocstring for ``from_delayed`` function  (:pr:`8576`) `Kirito1397`_
- Handle ``plot_width`` / ``plot_height`` deprecations (:pr:`8544`) `Bryan Van de Ven`_
- Remove unnecessary ``pyyaml`` ``importorskip`` (:pr:`8562`) `James Bourbeau`_
- Specify scheduler in DataFrame ``assert_eq`` (:pr:`8559`) `Gabe Joseph`_


.. _v2022.01.0:

2022.01.0
---------

Released on January 14, 2022

New Features
^^^^^^^^^^^^
- Add ``groupby.shift`` method (:pr:`8522`) `kori73`_
- Add ``DataFrame.nunique`` (:pr:`8479`) `Sarah Charlotte Johnson`_
- Add ``da.ndim`` to match ``np.ndim`` (:pr:`8502`) `Julia Signell`_

Enhancements
^^^^^^^^^^^^
- Only show ``percentile`` ``interpolation=`` keyword warning if NumPy version >= 1.22 (:pr:`8564`) `Julia Signell`_
- Raise ``PerformanceWarning`` when ``limit`` and ``"array.slicing.split-large-chunks"`` are ``None`` (:pr:`8511`) `Julia Signell`_
- Define ``normalize_seq`` function at import time (:pr:`8521`) `Illviljan`_
- Ensure that divisions are alway tuples (:pr:`8393`) `Charles Blackmon-Luca`_
- Allow a callable scheduler for ``bag.groupby`` (:pr:`8492`) `Julia Signell`_
- Save Zarr arrays with dask-on-ray scheduler (:pr:`8472`) `TnTo`_
- Make byte blocks more even in ``read_bytes`` (:pr:`8459`) `Martin Durant`_
- Improved the efficiency of ``matmul()`` by completely removing concatenation (:pr:`8423`) `ParticularMiner`_
- Limit max chunk size when reshaping dask arrays (:pr:`8124`) `Genevieve Buckley`_
- Changes for fastparquet superthrift (:pr:`8470`) `Martin Durant`_

Bug Fixes
^^^^^^^^^
- Fix boolean indices in array assignment (:pr:`8538`) `David Hassell`_
- Detect default ``dtype`` on array-likes (:pr:`8501`) `aeisenbarth`_
- Fix ``optimize_blockwise`` bug for duplicate dependency names (:pr:`8542`) `Richard (Rick) Zamora`_
- Update warnings for ``DataFrame.GroupBy.apply`` and transform (:pr:`8507`) `Sarah Charlotte Johnson`_
- Track HLG layer name in ``Delayed`` (:pr:`8452`) `Gabe Joseph`_
- Fix single item ``nanmin`` and ``nanmax`` reductions (:pr:`8484`) `Julia Signell`_
- Make ``read_csv`` with ``comment`` ``kwarg`` work even if there is a comment in the header (:pr:`8433`) `Julia Signell`_

Deprecations
^^^^^^^^^^^^
- Replace ``interpolation`` with ``method`` and ``method`` with ``internal_method`` (:pr:`8525`) `Julia Signell`_
- Remove daily stock demo utility (:pr:`8477`) `James Bourbeau`_

Documentation
^^^^^^^^^^^^^
- Add a join example in docs that be run with copy/paste (:pr:`8520`) `kori73`_
- Mention dashboard link in config (:pr:`8510`) `Ray Bell`_
- Fix changelog section hyperlinks (:pr:`8534`) `Aneesh Nema`_
- Hyphenate "single-machine scheduler" for consistency (:pr:`8519`) `Deepyaman Datta`_
- Normalize whitespace in doctests in ``slicing.py`` (:pr:`8512`) `Maren Westermann`_
- Best practices storage line typo (:pr:`8529`) `Michael Delgado`_
- Update figures (:pr:`8401`) `Sarah Charlotte Johnson`_
- Remove ``pyarrow``-only reference from ``split_row_groups`` in ``read_parquet`` docstring (:pr:`8490`) `Naty Clementi`_

Maintenance
^^^^^^^^^^^
- Remove obsolete ``LocalFileSystem`` tests that fail for ``fsspec>=2022.1.0`` (:pr:`8565`) `Richard (Rick) Zamora`_
- Tweak: "RuntimeWarning: invalid value encountered in reciprocal" (:pr:`8561`) `Guido Imperiale`_
- Fix ``skipna=None`` for ``DataFrame.sem`` (:pr:`8556`) `Julia Signell`_
- Fix ``PANDAS_GT_140`` (:pr:`8552`) `Julia Signell`_
- Collections with HLG must always implement ``__dask_layers__`` (:pr:`8548`) `Guido Imperiale`_
- Work around race condition in ``import llvmlite`` (:pr:`8550`) `Guido Imperiale`_
- Set a minimum version for ``pyyaml`` (:pr:`8545`) `Gaurav Sheni`_
- Adding ``nodefaults`` to environments to fix ``tiledb`` + mac issue (:pr:`8505`) `Julia Signell`_
- Set ceiling for ``setuptools`` (:pr:`8509`) `Julia Signell`_
- Add workflow / recipe to generate Dask nightlies (:pr:`8469`) `Charles Blackmon-Luca`_
- Bump gpuCI ``CUDA_VER`` to 11.5 (:pr:`8489`) `Charles Blackmon-Luca`_


.. _v2021.12.0:

2021.12.0
---------

Released on December 10, 2021

New Features
^^^^^^^^^^^^
- Add ``Series`` and ``Index`` ``is_monotonic*`` methods (:pr:`8304`) `Daniel Mesejo-León`_

Enhancements
^^^^^^^^^^^^
- Blockwise ``map_partitions`` with ``partition_info`` (:pr:`8310`) `Gabe Joseph`_
- Better error message for length of array with unknown chunk sizes (:pr:`8436`) `Doug Davis`_
- Use ``by`` instead of ``index`` internally on the Groupby class (:pr:`8441`) `Julia Signell`_
- Allow custom sort functions for ``sort_values`` (:pr:`8345`) `Charles Blackmon-Luca`_
- Add warning to ``read_parquet`` when statistics and partitions are misaligned (:pr:`8416`) `Richard (Rick) Zamora`_
- Support ``where`` argument in ufuncs (:pr:`8253`) `mihir`_
- Make visualize more consistent with compute (:pr:`8328`) `JSKenyon`_

Bug Fixes
^^^^^^^^^
- Fix ``map_blocks`` not using own arguments in ``name`` generation (:pr:`8462`) `David Hoese`_
- Fix for index error with reading empty parquet file (:pr:`8410`) `Sarah Charlotte Johnson`_
- Fix nullable-dtype error when writing partitioned parquet data (:pr:`8400`) `Richard (Rick) Zamora`_
- Fix CSV header bug (:pr:`8413`) `Richard (Rick) Zamora`_
- Fix empty chunk causes exception in ``nanmin``/``nanmax`` (:pr:`8375`) `Boaz Mohar`_

Deprecations
^^^^^^^^^^^^
- Deprecate ``token`` keyword argument to ``map_blocks`` (:pr:`8464`) `James Bourbeau`_
- Deprecation warning for default value of boundary kwarg in ``map_overlap`` (:pr:`8397`) `Genevieve Buckley`_

Documentation
^^^^^^^^^^^^^
- Clarify ``block_info`` documentation (:pr:`8425`) `Genevieve Buckley`_
- Output from alt text sprint (:pr:`8456`) `Sarah Charlotte Johnson`_
- Update talks and presentations (:pr:`8370`) `Naty Clementi`_
- Update Anaconda link in "Paid support" section of docs (:pr:`8427`) `Martin Durant`_
- Fixed broken ``dask-gateway`` link in ``ecosystem.rst`` (:pr:`8424`) `ofirr`_
- Fix CuPy doctest error (:pr:`8412`) `Genevieve Buckley`_

Maintenance
^^^^^^^^^^^
- Bump Bokeh min version to 2.1.1 (:pr:`8431`) `Bryan Van de Ven`_
- Fix following ``fsspec=2021.11.1`` release (:pr:`8428`) `Martin Durant`_
- Add ``dask/ml.py`` to pytest exclude list (:pr:`8414`) `Genevieve Buckley`_
- Update gpuCI ``RAPIDS_VER`` to ``22.02`` (:pr:`8394`)
- Unpin ``graphviz``  and improve package management in environment-3.7 (:pr:`8411`) `Julia Signell`_


.. _v2021.11.2:

2021.11.2
---------

Released on November 19, 2021

- Only run gpuCI bump script daily (:pr:`8404`) `Charles Blackmon-Luca`_
- Actually ignore index when asked in ``assert_eq`` (:pr:`8396`) `Gabe Joseph`_
- Ensure single-partition join ``divisions`` is ``tuple`` (:pr:`8389`) `Charles Blackmon-Luca`_
- Try to make divisions behavior clearer (:pr:`8379`) `Julia Signell`_
- Fix typo in ``set_index`` ``partition_size`` parameter description (:pr:`8384`) `FredericOdermatt`_
- Use ``blockwise`` in ``single_partition_join`` (:pr:`8341`) `Gabe Joseph`_
- Use more explicit keyword arguments (:pr:`8354`) `Boaz Mohar`_
- Fix ``.loc`` of DataFrame with nullable boolean ``dtype`` (:pr:`8368`) `Marco Rossi`_
- Parameterize shuffle implementation in tests (:pr:`8250`) `Ian Rose`_
- Remove some doc build warnings (:pr:`8369`) `Boaz Mohar`_
- Include properties in array API docs (:pr:`8356`) `Julia Signell`_
- Fix Zarr for upstream (:pr:`8367`) `Julia Signell`_
- Pin ``graphviz`` to avoid issue with windows and Python 3.7 (:pr:`8365`) `Julia Signell`_
- Import ``graphviz.Diagraph`` from top of module, not from ``dot`` (:pr:`8363`) `Julia Signell`_


.. _v2021.11.1:

2021.11.1
---------

Released on November 8, 2021

Patch release to update ``distributed`` dependency to version ``2021.11.1``.


.. _v2021.11.0:

2021.11.0
---------

Released on November 5, 2021

- Fx ``required_extension`` behavior in ``read_parquet`` (:pr:`8351`) `Richard (Rick) Zamora`_
- Add ``align_dataframes`` to ``map_partitions`` to broadcast a dataframe passed as an arg (:pr:`6628`) `Julia Signell`_
- Better handling for arrays/series of keys in ``dask.dataframe.loc`` (:pr:`8254`) `Julia Signell`_
- Point users to Discourse (:pr:`8332`) `Ian Rose`_
- Add ``name_function`` option to ``to_parquet`` (:pr:`7682`) `Matthew Powers`_
- Get rid of ``environment-latest.yml`` and update to Python 3.9 (:pr:`8275`) `Julia Signell`_
- Require newer ``s3fs`` in CI (:pr:`8336`) `James Bourbeau`_
- Groupby Rolling (:pr:`8176`) `Julia Signell`_
- Add more ordering diagnostics to ``dask.visualize`` (:pr:`7992`) `Erik Welch`_
- Use ``HighLevelGraph`` optimizations for ``delayed`` (:pr:`8316`) `Ian Rose`_
- ``demo_tuples`` produces malformed ``HighLevelGraph`` (:pr:`8325`) `Guido Imperiale`_
- Dask calendar should show events in local time (:pr:`8312`) `Genevieve Buckley`_
- Fix flaky ``test_interrupt`` (:pr:`8314`) `Guido Imperiale`_
- Deprecate ``AxisError`` (:pr:`8305`) `Guido Imperiale`_
- Fix name of cuDF in extension documentation. (:pr:`8311`) `Vyas Ramasubramani`_
- Add single eq operator (=) to parquet filters  (:pr:`8300`) `Ayush Dattagupta`_
- Improve support for Spark output in ``read_parquet`` (:pr:`8274`) `Richard (Rick) Zamora`_
- Add ``dask.ml`` module (:pr:`6384`) `Matthew Rocklin`_
- CI fixups (:pr:`8298`) `James Bourbeau`_
- Make slice errors match NumPy (:pr:`8248`) `Julia Signell`_
- Fix API docs misrendering with new sphinx theme (:pr:`8296`) `Julia Signell`_
- Replace ``block`` property with ``blockview`` for array-like operations on blocks (:pr:`8242`) `Davis Bennett`_
- Deprecate ``file_path`` and make it possible to save  from within a notebook (:pr:`8283`) `Julia Signell`_


.. _v2021.10.0:

2021.10.0
---------

Released on October 22, 2021

- ``da.store`` to create well-formed ``HighLevelGraph`` (:pr:`8261`) `Guido Imperiale`_
- CI: force nightly ``pyarrow`` in the upstream build (:pr:`8281`) `Joris Van den Bossche`_
- Remove ``chest`` (:pr:`8279`) `James Bourbeau`_
- Skip doctests if optional dependencies are not installed (:pr:`8258`) `Genevieve Buckley`_
- Update ``tmpdir`` and ``tmpfile`` context manager docstrings (:pr:`8270`) `Daniel Mesejo-León`_
- Unregister callbacks in doctests (:pr:`8276`) `James Bourbeau`_
- Fix typo in docs (:pr:`8277`) `JoranDox`_
- Stale label GitHub action (:pr:`8244`) `Genevieve Buckley`_
- Client-shutdown method appears twice (:pr:`8273`) `German Shiklov`_
- Add pre-commit to test requirements (:pr:`8257`) `Genevieve Buckley`_
- Refactor ``read_metadata`` in ``fastparquet`` engine (:pr:`8092`) `Richard (Rick) Zamora`_
- Support ``Path`` objects in ``from_zarr`` (:pr:`8266`) `Samuel Gaist`_
- Make nested redirects work (:pr:`8272`) `Julia Signell`_
- Set ``memory_usage`` to ``True`` if ``verbose`` is ``True`` in info (:pr:`8222`) `Kinshuk Dua`_
- Remove individual API doc pages from sphinx toctree (:pr:`8238`) `James Bourbeau`_
- Ignore whitespace in gufunc ``signature`` (:pr:`8267`) `James Bourbeau`_
- Add workflow to update gpuCI (:pr:`8215`) `Charles Blackmon-Luca`_
- ``DataFrame.head`` shouldn't warn when there's one partition (:pr:`8091`) `Pankaj Patil`_
- Ignore arrow doctests if ``pyarrow`` not installed (:pr:`8256`) `Genevieve Buckley`_
- Fix ``debugging.html`` redirect (:pr:`8251`) `James Bourbeau`_
- Fix null sorting for single partition dataframes (:pr:`8225`) `Charles Blackmon-Luca`_
- Fix ``setup.html`` redirect (:pr:`8249`) `Florian Jetter`_
- Run ``pyupgrade`` in CI (:pr:`8246`) `Guido Imperiale`_
- Fix label typo in upstream CI build (:pr:`8237`) `James Bourbeau`_
- Add support for "dependent" columns in DataFrame.assign (:pr:`8086`) `Suriya Senthilkumar`_
- add NumPy array of Dask keys to ``Array`` (:pr:`7922`) `Davis Bennett`_
- Remove unnecessary ``dask.multiprocessing`` import in docs (:pr:`8240`) `Ray Bell`_
- Adjust retrieving ``_max_workers`` from ``Executor`` (:pr:`8228`) `John A Kirkham`_
- Update function signatures in ``delayed`` best practices docs (:pr:`8231`) `Vũ Trung Đức`_
- Docs reoganization (:pr:`7984`) `Julia Signell`_
- Fix ``df.quantile`` on all missing data (:pr:`8129`) `Julia Signell`_
- Add ``tokenize.ensure-deterministic`` config option (:pr:`7413`) `Hristo Georgiev`_
- Use ``inclusive`` rather than ``closed`` with ``pandas>=1.4.0`` and ``pd.date_range`` (:pr:`8213`) `Julia Signell`_
- Add ``dask-gateway``, Coiled, and Saturn-Cloud to list of Dask setup tools (:pr:`7814`) `Kristopher Overholt`_
- Ensure existing futures get passed as deps when serializing ``HighLevelGraph`` layers (:pr:`8199`) `Jim Crist-Harif`_
- Make sure that the divisions of the single partition merge is left (:pr:`8162`) `Julia Signell`_
- Refactor ``read_metadata`` in ``pyarrow`` parquet engines (:pr:`8072`) `Richard (Rick) Zamora`_
- Support negative ``drop_axis`` in ``map_blocks`` and ``map_overlap`` (:pr:`8192`) `Gregory R. Lee`_
- Fix upstream tests (:pr:`8205`) `Julia Signell`_
- Add support for scalar item assignment by Series (:pr:`8195`) `Charles Blackmon-Luca`_
- Add some basic examples to doc strings on ``dask.bag`` ``all``, ``any``, ``count`` methods (:pr:`7630`) `Nathan Danielsen`_
- Don't have upstream report depend on commit message (:pr:`8202`) `James Bourbeau`_
- Ensure upstream CI cron job runs (:pr:`8200`) `James Bourbeau`_
- Use ``pytest.param`` to properly label param-specific GPU tests (:pr:`8197`) `Charles Blackmon-Luca`_
- Add ``test_set_index`` to tests ran on gpuCI (:pr:`8198`) `Charles Blackmon-Luca`_
- Suppress ``tmpfile`` OSError (:pr:`8191`) `James Bourbeau`_
- Use ``s.isna`` instead of ``pd.isna(s)`` in ``set_partitions_pre`` (fix cudf CI) (:pr:`8193`) `Charles Blackmon-Luca`_
- Open an issue for ``test-upstream`` failures (:pr:`8067`) `Wallace Reis`_
- Fix ``to_parquet`` bug in call to ``pyarrow.parquet.read_metadata`` (:pr:`8186`) `Richard (Rick) Zamora`_
- Add handling for null values in ``sort_values`` (:pr:`8167`) `Charles Blackmon-Luca`_
- Bump ``RAPIDS_VER`` for gpuCI (:pr:`8184`) `Charles Blackmon-Luca`_
- Dispatch walks MRO for lazily registered handlers (:pr:`8185`) `Jim Crist-Harif`_
- Configure SSHCluster instructions (:pr:`8181`) `Ray Bell`_
- Preserve ``HighLevelGraphs`` in ``DataFrame.from_delayed`` (:pr:`8174`) `Gabe Joseph`_
- Deprecate ``inplace`` argument for Dask series renaming (:pr:`8136`) `Marcel Coetzee`_
- Fix rolling for compatibility with ``pandas > 1.3.0`` (:pr:`8150`) `Julia Signell`_
- Raise error when ``setitem`` on unknown chunks (:pr:`8166`) `Julia Signell`_
- Include divisions when doing ``Index.to_series`` (:pr:`8165`) `Julia Signell`_


.. _v2021.09.1:

2021.09.1
---------

Released on September 21, 2021

- Fix ``groupby`` for future pandas (:pr:`8151`) `Julia Signell`_
- Remove warning filters in tests that are no longer needed (:pr:`8155`) `Julia Signell`_
- Add link to diagnostic visualize function in local diagnostic docs (:pr:`8157`) `David Hoese`_
- Add ``datetime_is_numeric`` to ``dataframe.describe`` (:pr:`7719`) `Julia Signell`_
- Remove references to ``pd.Int64Index`` in anticipation of deprecation (:pr:`8144`) `Julia Signell`_
- Use ``loc`` if needed for series ``__get_item__`` (:pr:`7953`) `Julia Signell`_
- Specifically ignore warnings on mean for empty slices (:pr:`8125`) `Julia Signell`_
- Skip ``groupby`` ``nunique`` test for pandas >= 1.3.3 (:pr:`8142`) `Julia Signell`_
- Implement ``ascending`` arg for ``sort_values`` (:pr:`8130`) `Charles Blackmon-Luca`_
- Replace ``operator.getitem``  (:pr:`8015`) `Naty Clementi`_
- Deprecate ``zero_broadcast_dimensions`` and ``homogeneous_deepmap`` (:pr:`8134`) `SnkSynthesis`_
- Add error if ``drop_index`` is negative (:pr:`8064`) `neel iyer`_
- Allow ``scheduler`` to be an ``Executor`` (:pr:`8112`) `John A Kirkham`_
- Handle ``asarray``/``asanyarray`` cases where ``like`` is a ``dask.Array`` (:pr:`8128`) `Peter Andreas Entschev`_
- Fix ``index_col`` duplication if ``index_col`` is type ``str`` (:pr:`7661`) `McToel`_
- Add ``dtype`` and ``order`` to ``asarray`` and ``asanyarray`` definitions (:pr:`8106`) `Julia Signell`_
- Deprecate ``dask.dataframe.Series.__contains__`` (:pr:`7914`) `Julia Signell`_
- Fix edge case with ``like``-arrays in ``_wrapped_qr`` (:pr:`8122`) `Peter Andreas Entschev`_
- Deprecate ``boundary_slice`` kwarg: ``kind`` for pandas compat (:pr:`8037`) `Julia Signell`_


.. _v2021.09.0:

2021.09.0
---------

Released on September 3, 2021

- Fewer open files (:pr:`7303`) `Julia Signell`_
- Add ``FileNotFound`` to expected http errors (:pr:`8109`) `Martin Durant`_
- Add ``DataFrame.sort_values`` to API docs (:pr:`8107`) `Benjamin Zaitlen`_
- Change to ``dask.order``: be more eager at times (:pr:`7929`) `Erik Welch`_
- Add pytest color to CI (:pr:`8090`) `James Bourbeau`_
- FIX: ``make_people`` works with ``processes`` scheduler (:pr:`8103`) `Dahn`_
- Adds ``deep`` param to Dataframe copy method and restrict it to ``False`` (:pr:`8068`) `João Paulo Lacerda`_
- Fix typo in configuration docs (:pr:`8104`) `Robert Hales`_
- Update formatting in ``DataFrame.query`` docstring (:pr:`8100`) `James Bourbeau`_
- Un-xfail ``sparse`` tests for 0.13.0 release (:pr:`8102`) `James Bourbeau`_
- Add axes property to DataFrame and Series (:pr:`8069`) `Jordan Jensen`_
- Add CuPy support in ``da.unique`` (values only) (:pr:`8021`) `Peter Andreas Entschev`_
- Unit tests for ``sparse.zeros_like`` (xfailed) (:pr:`8093`) `Guido Imperiale`_
- Add explicit ``like`` kwarg support to array creation functions (:pr:`8054`) `Peter Andreas Entschev`_
- Separate Array and DataFrame mindeps builds (:pr:`8079`) `James Bourbeau`_
- Fork out ``percentile_dispatch`` to ``dask.array`` (:pr:`8083`) `GALI PREM SAGAR`_
- Ensure ``filepath`` exists in ``to_parquet`` (:pr:`8057`) `James Bourbeau`_
- Update scheduler plugin usage in ``test_scheduler_highlevel_graph_unpack_import`` (:pr:`8080`) `James Bourbeau`_
- Add ``DataFrame.shuffle`` to API docs (:pr:`8076`) `Martin Fleischmann`_
- Order requirements alphabetically (:pr:`8073`) `John A Kirkham`_


.. _v2021.08.1:

2021.08.1
---------

Released on August 20, 2021

- Add ``ignore_metadata_file`` option to ``read_parquet`` (``pyarrow-dataset`` and ``fastparquet`` support only) (:pr:`8034`) `Richard (Rick) Zamora`_
- Add reference to ``pytest-xdist`` in dev docs (:pr:`8066`) `Julia Signell`_
- Include ``tz`` in meta from ``to_datetime`` (:pr:`8000`) `Julia Signell`_
- CI Infra Docs (:pr:`7985`) `Benjamin Zaitlen`_
- Include invalid DataFrame key in ``assert_eq`` check (:pr:`8061`) `James Bourbeau`_
- Use ``__class__`` when creating DataFrames (:pr:`8053`) `Mads R. B. Kristensen`_
- Use development version of ``distributed`` in gpuCI build (:pr:`7976`) `James Bourbeau`_
- Ignore whitespace when gufunc ``signature`` (:pr:`8049`) `James Bourbeau`_
- Move pandas import and percentile dispatch refactor (:pr:`8055`) `GALI PREM SAGAR`_
- Add colors to represent high level layer types (:pr:`7974`) `Freyam Mehta`_
- Upstream instance fix (:pr:`8060`) `Jacob Tomlinson`_
- Add ``dask.widgets`` and migrate HTML reprs to ``jinja2`` (:pr:`8019`) `Jacob Tomlinson`_
- Remove ``wrap_func_like_safe``, not required with NumPy >= 1.17 (:pr:`8052`) `Peter Andreas Entschev`_
- Fix threaded scheduler memory backpressure regression (:pr:`8040`) `David Hoese`_
- Add percentile dispatch (:pr:`8029`) `GALI PREM SAGAR`_
- Use a publicly documented attribute ``obj`` in ``groupby`` rather than private ``_selected_obj`` (:pr:`8038`) `GALI PREM SAGAR`_
- Specify module to ``import rechunk`` from (:pr:`8039`) `Illviljan`_
- Use ``dict`` to store data for {nan,}arg{min,max} in certain cases (:pr:`8014`) `Peter Andreas Entschev`_
- Fix ``blocksize`` description formatting in ``read_pandas`` (:pr:`8047`) `Louis Maddox`_
- Fix "point" -> "pointers" typo in docs (:pr:`8043`) `David Chudzicki`_


.. _v2021.08.0:

2021.08.0
---------

Released on August 13, 2021

- Fix ``to_orc`` delayed compute behavior (:pr:`8035`) `Richard (Rick) Zamora`_
- Don't convert to low-level task graph in ``compute_as_if_collection`` (:pr:`7969`) `James Bourbeau`_
- Fix multifile read for hdf (:pr:`8033`) `Julia Signell`_
- Resolve warning in ``distributed`` tests (:pr:`8025`) `James Bourbeau`_
- Update ``to_orc`` collection name (:pr:`8024`) `James Bourbeau`_
- Resolve ``skipfooter`` problem (:pr:`7855`) `Ross`_
- Raise ``NotImplementedError`` for non-indexable arg passed to ``to_datetime`` (:pr:`7989`) `Doug Davis`_
- Ensure we error on warnings from ``distributed`` (:pr:`8002`) `James Bourbeau`_
- Added ``dict`` format in ``to_bag`` accessories of DataFrame (:pr:`7932`) `gurunath`_
- Delayed docs indirect dependencies (:pr:`8016`) `aa1371`_
- Add tooltips to graphviz high-level graphs (:pr:`7973`) `Freyam Mehta`_
- Close 2021 User Survey (:pr:`8007`) `Julia Signell`_
- Reorganize CuPy tests into multiple files (:pr:`8013`) `Peter Andreas Entschev`_
- Refactor and Expand Dask-Dataframe ORC API (:pr:`7756`) `Richard (Rick) Zamora`_
- Don't enforce columns if ``enforce=False`` (:pr:`7916`) `Julia Signell`_
- Fix ``map_overlap`` trimming behavior when ``drop_axis`` is not ``None`` (:pr:`7894`) `Gregory R. Lee`_
- Mark gpuCI CuPy test as flaky (:pr:`7994`) `Peter Andreas Entschev`_
- Avoid using ``Delayed`` in ``to_csv`` and ``to_parquet`` (:pr:`7968`) `Matthew Rocklin`_
- Removed redundant ``check_dtypes`` (:pr:`7952`) `gurunath`_
- Use ``pytest.warns`` instead of raises for checking parquet engine deprecation (:pr:`7993`) `Joris Van den Bossche`_
- Bump ``RAPIDS_VER`` in gpuCI to 21.10 (:pr:`7991`) `Charles Blackmon-Luca`_
- Add back ``pyarrow-legacy`` test coverage for ``pyarrow>=5`` (:pr:`7988`) `Richard (Rick) Zamora`_
- Allow ``pyarrow>=5`` in ``to_parquet`` and ``read_parquet`` (:pr:`7967`) `Richard (Rick) Zamora`_
- Skip CuPy tests requiring NEP-35 when NumPy < 1.20 is available (:pr:`7982`) `Peter Andreas Entschev`_
- Add ``tail`` and ``head`` to ``SeriesGroupby`` (:pr:`7935`) `Daniel Mesejo-León`_
- Update Zoom link for monthly meeting (:pr:`7979`) `James Bourbeau`_
- Add gpuCI build script (:pr:`7966`) `Charles Blackmon-Luca`_
- Deprecate ``daily_stock`` utility (:pr:`7949`) `James Bourbeau`_
- Add ``distributed.nanny`` to configuration reference docs (:pr:`7955`) `James Bourbeau`_
- Require NumPy 1.18+ & Pandas 1.0+ (:pr:`7939`) `John A Kirkham`_


.. _v2021.07.2:

2021.07.2
---------

Released on July 30, 2021

.. note::

  This is the last release with support for NumPy 1.17 and pandas 0.25.
  Beginning with the next release, NumPy 1.18 and pandas 1.0 will be the minimum
  supported versions.

- Add ``dask.array`` SVG to the HTML Repr (:pr:`7886`) `Freyam Mehta`_
- Avoid use of ``Delayed`` in ``to_parquet`` (:pr:`7958`) `Matthew Rocklin`_
- Temporarily pin ``pyarrow<5`` in CI (:pr:`7960`) `James Bourbeau`_
- Add deprecation warning for top-level ``ucx`` and ``rmm`` config values (:pr:`7956`) `James Bourbeau`_
- Remove skips from doctests (4 of 6) (:pr:`7865`) `Zhengnan Zhao`_
- Remove skips from doctests (5 of 6) (:pr:`7864`) `Zhengnan Zhao`_
- Adds missing prepend/append functionality to ``da.diff`` (:pr:`7946`) `Peter Andreas Entschev`_
- Change graphviz font family to sans (:pr:`7931`) `Freyam Mehta`_
- Fix read-csv name - when path is different, use different name for task (:pr:`7942`) `Julia Signell`_
- Update configuration reference for ``ucx`` and ``rmm`` changes (:pr:`7943`) `James Bourbeau`_
- Add meta support to ``__setitem__`` (:pr:`7940`) `Peter Andreas Entschev`_
- NEP-35 support for ``slice_with_int_dask_array`` (:pr:`7927`) `Peter Andreas Entschev`_
- Unpin fastparquet in CI (:pr:`7928`) `James Bourbeau`_
- Remove skips from doctests (3 of 6) (:pr:`7872`) `Zhengnan Zhao`_


.. _v2021.07.1:

2021.07.1
---------

Released on July 23, 2021

- Make array ``assert_eq`` check dtype (:pr:`7903`) `Julia Signell`_
- Remove skips from doctests (6 of 6) (:pr:`7863`) `Zhengnan Zhao`_
- Remove experimental feature warning from actors docs (:pr:`7925`) `Matthew Rocklin`_
- Remove skips from doctests (2 of 6) (:pr:`7873`) `Zhengnan Zhao`_
- Separate out Array and Bag API (:pr:`7917`) `Julia Signell`_
- Implement lazy ``Array.__iter__`` (:pr:`7905`) `Julia Signell`_
- Clean up places where we inadvertently iterate over arrays (:pr:`7913`) `Julia Signell`_
- Add ``numeric_only`` kwarg to DataFrame reductions (:pr:`7831`) `Julia Signell`_
- Add pytest marker for GPU tests (:pr:`7876`) `Charles Blackmon-Luca`_
- Add support for ``histogram2d`` in ``dask.array`` (:pr:`7827`) `Doug Davis`_
- Remove skips from doctests (1 of 6) (:pr:`7874`) `Zhengnan Zhao`_
- Add node size scaling to the Graphviz output for the high level graphs (:pr:`7869`) `Freyam Mehta`_
- Update old Bokeh links (:pr:`7915`) `Bryan Van de Ven`_
- Temporarily pin ``fastparquet`` in CI (:pr:`7907`) `James Bourbeau`_
- Add ``dask.array`` import to progress bar docs (:pr:`7910`) `Fabian Gebhart`_
- Use separate files for each DataFrame API function and method (:pr:`7890`) `Julia Signell`_
- Fix ``pyarrow-dataset`` ordering bug (:pr:`7902`) `Richard (Rick) Zamora`_
- Generalize unique aggregate (:pr:`7892`) `GALI PREM SAGAR`_
- Raise ``NotImplementedError`` when using ``pd.Grouper`` (:pr:`7857`) `Ruben van de Geer`_
- Add ``aggregate_files`` argument to enable multi-file partitions in ``read_parquet`` (:pr:`7557`) `Richard (Rick) Zamora`_
- Un-``xfail`` ``test_daily_stock`` (:pr:`7895`) `James Bourbeau`_
- Update access configuration docs (:pr:`7837`) `Naty Clementi`_
- Use packaging for version comparisons (:pr:`7820`) `Elliott Sales de Andrade`_
- Handle infinite loops in ``merge_asof`` (:pr:`7842`) `gerrymanoim`_


.. _v2021.07.0:

2021.07.0
---------

Released on July 9, 2021

- Include ``fastparquet`` in upstream CI build (:pr:`7884`) `James Bourbeau`_
- Blockwise: handle non-string constant dependencies  (:pr:`7849`) `Mads R. B. Kristensen`_
- ``fastparquet`` now supports new time types, including ns precision (:pr:`7880`) `Martin Durant`_
- Avoid ``ParquetDataset`` API when appending in ``ArrowDatasetEngine`` (:pr:`7544`) `Richard (Rick) Zamora`_
- Add retry logic to ``test_shuffle_priority`` (:pr:`7879`) `Richard (Rick) Zamora`_
- Use strict channel priority in CI (:pr:`7878`) `James Bourbeau`_
- Support nested ``dask.distributed`` imports (:pr:`7866`) `Matthew Rocklin`_
- Should check module name only, not the entire directory filepath (:pr:`7856`) `Genevieve Buckley`_
- Updates due to https://github.com/dask/fastparquet/pull/623 (:pr:`7875`) `Martin Durant`_
- ``da.eye`` fix for ``chunks=-1`` (:pr:`7854`) `Naty Clementi`_
- Temporarily xfail ``test_daily_stock`` (:pr:`7858`) `James Bourbeau`_
- Set priority annotations in ``SimpleShuffleLayer`` (:pr:`7846`) `Richard (Rick) Zamora`_
- Blockwise: stringify constant key inputs (:pr:`7838`) `Mads R. B. Kristensen`_
- Allow mixing dask and numpy arrays in ``@guvectorize`` (:pr:`6863`) `Julia Signell`_
- Don't sample dict result of a shuffle group when calculating its size (:pr:`7834`) `Florian Jetter`_
- Fix scipy tests (:pr:`7841`) `Julia Signell`_
- Deterministically tokenize ``datetime.date`` (:pr:`7836`) `James Bourbeau`_
- Add ``sample_rows`` to ``read_csv``-like (:pr:`7825`) `Martin Durant`_
- Fix typo in ``config.deserialize`` docstring (:pr:`7830`) `Geoffrey Lentner`_
- Remove warning filter in ``test_dataframe_picklable`` (:pr:`7822`) `James Bourbeau`_
- Improvements to ``histogramdd`` (for handling inputs that are sequences-of-arrays). (:pr:`7634`) `Doug Davis`_
- Make ``PY_VERSION`` private (:pr:`7824`) `James Bourbeau`_


.. _v2021.06.2:

2021.06.2
---------

Released on June 22, 2021

- ``layers.py`` compare ``parts_out`` with ``set(self.parts_out)`` (:pr:`7787`) `Genevieve Buckley`_
- Make ``check_meta`` understand pandas dtypes better (:pr:`7813`) `Julia Signell`_
- Remove "Educational Resources" doc page (:pr:`7818`) `James Bourbeau`_


.. _v2021.06.1:

2021.06.1
---------

Released on June 18, 2021

- Replace funding page with 'Supported By' section on dask.org (:pr:`7817`) `James Bourbeau`_
- Add initial deprecation utilities (:pr:`7810`) `James Bourbeau`_
- Enforce dtype conservation in ufuncs that explicitly use ``dtype=`` (:pr:`7808`) `Doug Davis`_
- Add Coiled to list of paid support organizations (:pr:`7811`) `Kristopher Overholt`_
- Small tweaks to the HTML repr for ``Layer`` & ``HighLevelGraph`` (:pr:`7812`) `Genevieve Buckley`_
- Add dark mode support to HLG HTML repr (:pr:`7809`) `Jacob Tomlinson`_
- Remove compatibility entries for old distributed (:pr:`7801`) `Elliott Sales de Andrade`_
- Implementation of HTML repr for ``HighLevelGraph`` layers (:pr:`7763`) `Genevieve Buckley`_
- Update default ``blockwise`` token to avoid DataFrame column name clash (:pr:`6546`) `James Bourbeau`_
- Use dispatch ``concat`` for ``merge_asof`` (:pr:`7806`) `Julia Signell`_
- Fix upstream freq tests (:pr:`7795`) `Julia Signell`_
- Use more context managers from the standard library (:pr:`7796`) `James Bourbeau`_
- Simplify skips in parquet tests (:pr:`7802`) `Elliott Sales de Andrade`_
- Remove check for outdated bokeh (:pr:`7804`) `Elliott Sales de Andrade`_
- More test coverage uploads (:pr:`7799`) `James Bourbeau`_
- Remove ``ImportError`` catching from ``dask/__init__.py`` (:pr:`7797`) `James Bourbeau`_
- Allow ``DataFrame.join()`` to take a list of DataFrames to merge with (:pr:`7578`) `Krishan Bhasin`_
- Fix maximum recursion depth exception in ``dask.array.linspace`` (:pr:`7667`) `Daniel Mesejo-León`_
- Fix docs links (:pr:`7794`) `Julia Signell`_
- Initial ``da.select()`` implementation and test (:pr:`7760`) `Gabriel Miretti`_
- Layers must implement ``get_output_keys`` method (:pr:`7790`) `Genevieve Buckley`_
- Don't include or expect ``freq`` in divisions (:pr:`7785`) `Julia Signell`_
- A ``HighLevelGraph`` abstract layer for ``map_overlap`` (:pr:`7595`) `Genevieve Buckley`_
- Always include kwarg name in ``drop`` (:pr:`7784`) `Julia Signell`_
- Only rechunk for median if needed (:pr:`7782`) `Julia Signell`_
- Add ``add_(prefix|suffix)`` to DataFrame and Series (:pr:`7745`) `tsuga`_
- Move ``read_hdf`` to ``Blockwise`` (:pr:`7625`) `Richard (Rick) Zamora`_
- Make ``Layer.get_output_keys`` officially an abstract method (:pr:`7775`) `Genevieve Buckley`_
- Non-dask-arrays and broadcasting in ``ravel_multi_index`` (:pr:`7594`) `Gabe Joseph`_
- Fix for paths ending with "/" in parquet overwrite (:pr:`7773`) `Martin Durant`_
- Fixing calling ``.visualize()`` with ``filename=None`` (:pr:`7740`) `Freyam Mehta`_
- Generate unique names for ``SubgraphCallable`` (:pr:`7637`) `Bruce Merry`_
- Pin ``fsspec`` to ``2021.5.0`` in CI (:pr:`7771`) `James Bourbeau`_
- Evaluate graph lazily if meta is provided in ``from_delayed`` (:pr:`7769`) `Florian Jetter`_
- Add ``meta`` support for ``DatetimeTZDtype`` (:pr:`7627`) `gerrymanoim`_
- Add dispatch label to automatic PR labeler (:pr:`7701`) `James Bourbeau`_
- Fix HDFS tests (:pr:`7752`) `Julia Signell`_


.. _v2021.06.0:

2021.06.0
---------

Released on June 4, 2021

- Remove abstract tokens from graph keys in ``rewrite_blockwise`` (:pr:`7721`) `Richard (Rick) Zamora`_
- Ensure correct column order in csv ``project_columns`` (:pr:`7761`) `Richard (Rick) Zamora`_
- Renamed inner loop variables to avoid duplication (:pr:`7741`) `Boaz Mohar`_
- Do not return delayed object from ``to_zarr`` (:pr:`7738`) `Chris Roat`
- Array: correct number of outputs in ``apply_gufunc`` (:pr:`7669`) `Gabe Joseph`_
- Rewrite ``da.fromfunction`` with ``da.blockwise`` (:pr:`7704`) `John A Kirkham`_
- Rename ``make_meta_util`` to ``make_meta`` (:pr:`7743`) `GALI PREM SAGAR`_
- Repartition before shuffle if the requested partitions are less than input partitions (:pr:`7715`) `Vibhu Jawa`_
- Blockwise: handle constant key inputs (:pr:`7734`) `Mads R. B. Kristensen`_
- Added raise to ``apply_gufunc`` (:pr:`7744`) `Boaz Mohar`_
- Show failing tests summary in CI (:pr:`7735`) `Genevieve Buckley`_
- ``sizeof`` sets in Python 3.9 (:pr:`7739`) `Mads R. B. Kristensen`_
- Warn if using pandas datetimelike string in ``dataframe.__getitem__`` (:pr:`7749`) `Julia Signell`_
- Highlight the ``client.dashboard_link`` (:pr:`7747`) `Genevieve Buckley`_
- Easier link for subscribing to the Google calendar (:pr:`7733`) `Genevieve Buckley`_
- Automatically show graph visualization in Jupyter notebooks (:pr:`7716`) `Genevieve Buckley`_
- Add ``autofunction`` for ``unify_chunks`` in API docs (:pr:`7730`) `James Bourbeau`_


.. _v2021.05.1:

2021.05.1
---------

Released on May 28, 2021

- Pandas compatibility (:pr:`7712`) `Julia Signell`_
- Fix ``optimize_dataframe_getitem`` bug (:pr:`7698`) `Richard (Rick) Zamora`_
- Update ``make_meta`` import in docs (:pr:`7713`) `Benjamin Zaitlen`_
- Implement ``da.searchsorted`` (:pr:`7696`) `Tom White`_
- Fix format string in error message (:pr:`7706`) `Jiaming Yuan`_
- Fix ``read_sql_table`` returning wrong result for single column loads (:pr:`7572`) `c-thiel`_
- Add slack join link in ``support.rst`` (:pr:`7679`) `Naty Clementi`_
- Remove unused alphabet variable (:pr:`7700`) `James Bourbeau`_
- Fix meta creation incase of ``object`` (:pr:`7586`) `GALI PREM SAGAR`_
- Add dispatch for ``union_categoricals`` (:pr:`7699`) `GALI PREM SAGAR`_
- Consolidate array ``Dispatch`` objects (:pr:`7505`) `James Bourbeau`_
- Move DataFrame ``dispatch.registers`` to their own file (:pr:`7503`) `Julia Signell`_
- Fix delayed with ``dataclasses`` where ``init=False`` (:pr:`7656`) `Julia Signell`_
- Allow a column to be named ``divisions`` (:pr:`7605`) `Julia Signell`_
- Stack nd array with unknown chunks (:pr:`7562`) `Chris Roat`_
- Promote the 2021 Dask User Survey (:pr:`7694`) `Genevieve Buckley`_
- Fix typo in ``DataFrame.set_index()`` (:pr:`7691`) `James Lamb`_
- Cleanup array API reference links (:pr:`7684`) `David Hoese`_
- Accept ``axis`` tuple for ``flip`` to be consistent with NumPy (:pr:`7675`) `Andrew Champion`_
- Bump ``pre-commit`` hook versions (:pr:`7676`) `James Bourbeau`_
- Cleanup ``to_zarr`` docstring (:pr:`7683`) `David Hoese`_
- Fix the docstring of ``read_orc`` (:pr:`7678`) `keewis`_
- Doc ``ipyparallel`` & ``mpi4py`` ``concurrent.futures`` (:pr:`7665`) `John A Kirkham`_
- Update tests to support CuPy 9 (:pr:`7671`) `Peter Andreas Entschev`_
- Fix some ``HighLevelGraph`` documentation inaccuracies (:pr:`7662`) `Mads R. B. Kristensen`_
- Fix spelling in Series ``getitem`` error message (:pr:`7659`) `Maisie Marshall`_


.. _v2021.05.0:

2021.05.0
---------

Released on May 14, 2021

- Remove deprecated ``kind`` kwarg to comply with pandas 1.3.0 (:pr:`7653`) `Julia Signell`_
- Fix bug in DataFrame column projection (:pr:`7645`) `Richard (Rick) Zamora`_
- Merge global annotations when packing (:pr:`7565`) `Mads R. B. Kristensen`_
- Avoid ``inplace=`` in pandas ``set_categories`` (:pr:`7633`) `James Bourbeau`_
- Change the active-fusion default to ``False`` for Dask-Dataframe (:pr:`7620`) `Richard (Rick) Zamora`_
- Array: remove extraneous code from ``RandomState`` (:pr:`7487`) `Gabe Joseph`_
- Implement ``str.concat`` when ``others=None`` (:pr:`7623`) `Daniel Mesejo-León`_
- Fix ``dask.dataframe`` in sandboxed environments (:pr:`7601`) `Noah D. Brenowitz`_
- Support for ``cupyx.scipy.linalg`` (:pr:`7563`) `Benjamin Zaitlen`_
- Move ``timeseries`` and daily-stock to ``Blockwise`` (:pr:`7615`) `Richard (Rick) Zamora`_
- Fix bugs in broadcast join (:pr:`7617`) `Richard (Rick) Zamora`_
- Use ``Blockwise`` for DataFrame IO (parquet, csv, and orc) (:pr:`7415`) `Richard (Rick) Zamora`_
- Adding chunk & type information to Dask ``HighLevelGraph`` s (:pr:`7309`) `Genevieve Buckley`_
- Add ``pyarrow`` sphinx ``intersphinx_mapping`` (:pr:`7612`) `Ray Bell`_
- Remove skip on test freq (:pr:`7608`) `Julia Signell`_
- Defaults in ``read_parquet`` parameters (:pr:`7567`) `Ray Bell`_
- Remove ``ignore_abc_warning`` (:pr:`7606`) `Julia Signell`_
- Harden DataFrame merge between column-selection and index (:pr:`7575`) `Richard (Rick) Zamora`_
- Get rid of ``ignore_abc`` decorator (:pr:`7604`) `Julia Signell`_
- Remove kwarg validation for bokeh (:pr:`7597`) `Julia Signell`_
- Add ``loky`` example (:pr:`7590`) `Naty Clementi`_
- Delayed: ``nout`` when arguments become tasks (:pr:`7593`) `Gabe Joseph`_
- Update distributed version in mindep CI build (:pr:`7602`) `James Bourbeau`_
- Support all or no overlap between partition columns and real columns (:pr:`7541`) `Richard (Rick) Zamora`_


.. _v2021.04.1:

2021.04.1
---------

Released on April 23, 2021

- Handle ``Blockwise`` HLG pack/unpack for ``concatenate=True`` (:pr:`7455`) `Richard (Rick) Zamora`_
- ``map_partitions``: use tokenized info as name of the ``SubgraphCallable`` (:pr:`7524`) `Mads R. B. Kristensen`_
- Using ``tmp_path`` and ``tmpdir`` to avoid temporary files and directories hanging in the repo (:pr:`7592`) `Naty Clementi`_
- Contributing to docs (development guide) (:pr:`7591`) `Naty Clementi`_
- Add more packages to Python 3.9 CI build (:pr:`7588`) `James Bourbeau`_
- Array: Fix NEP-18 dispatching in finalize (:pr:`7508`) `Gabe Joseph`_
- Misc fixes for ``numpydoc`` (:pr:`7569`) `Matthias Bussonnier`_
- Avoid pandas ``level=`` keyword deprecation (:pr:`7577`) `James Bourbeau`_
- Map e.g. ``.repartition(freq="M")`` to ``.repartition(freq="MS")`` (:pr:`7504`) `Ruben van de Geer`_
- Remove hash seeding in parallel CI runs (:pr:`7128`) `Elliott Sales de Andrade`_
- Add defaults in parameters in ``to_parquet`` (:pr:`7564`) `Ray Bell`_
- Simplify transpose axes cleanup (:pr:`7561`) `Julia Signell`_
- Make ``ValueError in len(index_names) > 1`` explicit it's using ``fastparquet`` (:pr:`7556`) `Ray Bell`_
- Fix ``dict``-column appending for ``pyarrow`` parquet engines (:pr:`7527`) `Richard (Rick) Zamora`_
- Add a documentation auto label (:pr:`7560`) `Doug Davis`_
- Add ``dask.delayed.Delayed`` to docs so it can be referenced by other sphinx docs (:pr:`7559`) `Doug Davis`_
- Fix upstream ``idxmaxmin`` for uneven ``split_every`` (:pr:`7538`) `Julia Signell`_
- Make ``normalize_token`` for pandas ``Series``/``DataFrame`` future proof (no direct block access) (:pr:`7318`) `Joris Van den Bossche`_
- Redesigned ``__setitem__`` implementation (:pr:`7393`) `David Hassell`_
- ``histogram``, ``histogramdd`` improvements (docs; return consistencies) (:pr:`7520`) `Doug Davis`_
- Force nightly ``pyarrow`` in the upstream build (:pr:`7530`) `Joris Van den Bossche`_
- Fix Configuration Reference (:pr:`7533`) `Benjamin Zaitlen`_
- Use ``.to_parquet`` on ``dask.dataframe`` in doc string (:pr:`7528`) `Ray Bell`_
- Avoid double ``msgpack`` serialization of HLGs (:pr:`7525`) `Mads R. B. Kristensen`_
- Encourage usage of ``yaml.safe_load()`` in configuration doc (:pr:`7529`) `Hristo Georgiev`_
- Fix ``reshape`` bug. Add relevant test. Fixes #7171. (:pr:`7523`) `JSKenyon`_
- Support ``custom_metadata=`` argument in ``to_parquet`` (:pr:`7359`) `Richard (Rick) Zamora`_
- Clean some documentation warnings (:pr:`7518`) `Daniel Mesejo-León`_
- Getting rid of more docs warnings (:pr:`7426`) `Julia Signell`_
- Added ``product`` (alias of ``prod``) (:pr:`7517`) `Freyam Mehta`_
- Fix upstream ``__array_ufunc__`` tests (:pr:`7494`) `Julia Signell`_
- Escape from ``map_overlap`` to ``map_blocks`` if depth is zero (:pr:`7481`) `Genevieve Buckley`_
- Add ``check_type`` to array ``assert_eq`` (:pr:`7491`) `Julia Signell`_


.. _v2021.04.0:

2021.04.0
---------

Released on April 2, 2021

- Adding support for multidimensional histograms with ``dask.array.histogramdd`` (:pr:`7387`) `Doug Davis`_
- Update docs on number of threads and workers in default ``LocalCluster`` (:pr:`7497`) `cameron16`_
- Add labels automatically when certain files are touched in a PR (:pr:`7506`) `Julia Signell`_
- Extract ``ignore_order`` from ``kwargs`` (:pr:`7500`) `GALI PREM SAGAR`_
- Only provide installation instructions when distributed is missing (:pr:`7498`) `Matthew Rocklin`_
- Start adding ``isort`` (:pr:`7370`) `Julia Signell`_
- Add ``ignore_order`` parameter in ``dd.concat`` (:pr:`7473`) `Daniel Mesejo-León`_
- Use powers-of-two when displaying RAM (:pr:`7484`) `Guido Imperiale`_
- Added License Classifier (:pr:`7485`) `Tom Augspurger`_
- Replace conda with mamba (:pr:`7227`) `Guido Imperiale`_
- Fix typo in array docs (:pr:`7478`) `James Lamb`_
- Use ``concurrent.futures`` in local scheduler (:pr:`6322`) `John A Kirkham`_


.. _v2021.03.1:

2021.03.1
---------

Released on March 26, 2021

- Add a dispatch for ``is_categorical_dtype`` to handle non-pandas objects (:pr:`7469`) `brandon-b-miller`_
- Use ``multiprocessing.Pool`` in ``test_read_text`` (:pr:`7472`) `John A Kirkham`_
- Add missing ``meta`` kwarg to gufunc class (:pr:`7423`) `Peter Andreas Entschev`_
- Example for memory-mapped Dask array (:pr:`7380`) `Dieter Weber`_
- Fix NumPy upstream failures ``xfail`` pandas and fastparquet failures (:pr:`7441`) `Julia Signell`_
- Fix bug in repartition with freq (:pr:`7357`) `Ruben van de Geer`_
- Fix ``__array_function__`` dispatching for ``tril``/``triu`` (:pr:`7457`) `Peter Andreas Entschev`_
- Use ``concurrent.futures.Executors`` in a few tests (:pr:`7429`) `John A Kirkham`_
- Require NumPy >=1.16 (:pr:`7383`) `Guido Imperiale`_
- Minor ``sort_values`` housekeeping (:pr:`7462`) `Ryan Williams`_
- Ensure natural sort order in parquet part paths (:pr:`7249`) `Ryan Williams`_
- Remove global env mutation upon running ``test_config.py`` (:pr:`7464`) `Hristo Georgiev`_
- Update NumPy intersphinx URL (:pr:`7460`) `Gabe Joseph`_
- Add ``rot90`` (:pr:`7440`) `Trevor Manz`_
- Update docs for required package for endpoint (:pr:`7454`) `Nick Vazquez`_
- Master -> main in ``slice_array`` docstring (:pr:`7453`) `Gabe Joseph`_
- Expand ``dask.utils.is_arraylike`` docstring (:pr:`7445`) `Doug Davis`_
- Simplify ``BlockwiseIODeps`` importing (:pr:`7420`) `Richard (Rick) Zamora`_
- Update layer annotation packing method (:pr:`7430`) `James Bourbeau`_
- Drop duplicate test in ``test_describe_empty`` (:pr:`7431`) `John A Kirkham`_
- Add ``Series.dot`` method to dataframe module (:pr:`7236`) `Madhu94`_
- Added df ``kurtosis``-method and testing (:pr:`7273`) `Jan Borchmann`_
- Avoid quadratic-time performance for HLG culling (:pr:`7403`) `Bruce Merry`_
- Temporarily skip problematic ``sparse`` test (:pr:`7421`) `James Bourbeau`_
- Update some CI workflow names (:pr:`7422`) `James Bourbeau`_
- Fix HDFS test (:pr:`7418`) `Julia Signell`_
- Make changelog subtitles match the hierarchy (:pr:`7419`) `Julia Signell`_
- Add support for normalize in ``value_counts`` (:pr:`7342`) `Julia Signell`_
- Avoid unnecessary imports for HLG Layer unpacking and materialization (:pr:`7381`) `Richard (Rick) Zamora`_
- Bincount fix slicing (:pr:`7391`) `Genevieve Buckley`_
- Add ``sliding_window_view`` (:pr:`7234`) `Deepak Cherian`_
- Fix typo in ``docs/source/develop.rst`` (:pr:`7414`) `Hristo Georgiev`_
- Switch documentation builds for PRs to readthedocs (:pr:`7397`) `James Bourbeau`_
- Adds ``sort_values`` to dask.DataFrame (:pr:`7286`) `gerrymanoim`_
- Pin ``sqlalchemy<1.4.0`` in CI (:pr:`7405`) `James Bourbeau`_
- Comment fixes (:pr:`7215`) `Ryan Williams`_
- Dead code removal / fixes (:pr:`7388`) `Ryan Williams`_
- Use single thread for ``pa.Table.from_pandas`` calls (:pr:`7347`) `Richard (Rick) Zamora`_
- Replace ``'container'`` with ``'image'`` (:pr:`7389`) `James Lamb`_
- DOC hyperlink repartition (:pr:`7394`) `Ray Bell`_
- Pass delimiter to ``fsspec`` in ``bag.read_text`` (:pr:`7349`) `Martin Durant`_
- Update ``read_hdf`` default mode to ``"r"`` (:pr:`7039`) `rs9w33`_
- Embed literals in ``SubgraphCallable`` when packing ``Blockwise`` (:pr:`7353`) `Mads R. B. Kristensen`_
- Update ``test_hdf.py`` to not reuse file handlers (:pr:`7044`) `rs9w33`_
- Require additional dependencies: cloudpickle, partd, fsspec, toolz (:pr:`7345`) `Julia Signell`_
- Prepare ``Blockwise`` + IO infrastructure (:pr:`7281`) `Richard (Rick) Zamora`_
- Remove duplicated imports from ``test_slicing.py`` (:pr:`7365`) `Hristo Georgiev`_
- Add test deps for pip development (:pr:`7360`) `Julia Signell`_
- Support int slicing for non-NumPy arrays (:pr:`7364`) `Peter Andreas Entschev`_
- Automatically cancel previous CI builds (:pr:`7348`) `James Bourbeau`_
- ``dask.array.asarray`` should handle case where ``xarray`` class is in top-level namespace (:pr:`7335`) `Tom White`_
- ``HighLevelGraph`` length without materializing layers (:pr:`7274`) `Gabe Joseph`_
- Drop support for Python 3.6 (:pr:`7006`) `James Bourbeau`_
- Fix fsspec usage in ``create_metadata_file`` (:pr:`7295`) `Richard (Rick) Zamora`_
- Change default branch from master to main (:pr:`7198`) `Julia Signell`_
- Add Xarray to CI software environment (:pr:`7338`) `James Bourbeau`_
- Update repartition argument name in error text (:pr:`7336`) `Eoin Shanaghy`_
- Run upstream tests based on commit message (:pr:`7329`) `James Bourbeau`_
- Use ``pytest.register_assert_rewrite`` on util modules (:pr:`7278`) `Bruce Merry`_
- Add example on using specific chunk sizes in ``from_array()`` (:pr:`7330`) `James Lamb`_
- Move NumPy skip into test (:pr:`7247`) `Julia Signell`_


.. _v2021.03.0:

2021.03.0
---------

Released on March 5, 2021

.. note::

    This is the first release with support for Python 3.9 and the
    last release with support for Python 3.6

- Bump minimum version of ``distributed`` (:pr:`7328`) `James Bourbeau`_
- Fix ``percentiles_summary`` with ``dask_cudf`` (:pr:`7325`) `Peter Andreas Entschev`_
- Temporarily revert recent ``Array.__setitem__`` updates (:pr:`7326`) `James Bourbeau`_
- ``Blockwise.clone`` (:pr:`7312`) `Guido Imperiale`_
- NEP-35 duck array update (:pr:`7321`) `James Bourbeau`_
- Don't allow setting ``.name`` for array (:pr:`7222`) `Julia Signell`_
- Use nearest interpolation for creating percentiles of integer input (:pr:`7305`) `Kyle Barron`_
- Test ``exp`` with CuPy arrays (:pr:`7322`) `John A Kirkham`_
- Check that computed chunks have right size and dtype (:pr:`7277`) `Bruce Merry`_
- ``pytest.mark.flaky`` (:pr:`7319`) `Guido Imperiale`_
- Contributing docs: add note to pull the latest git tags before pip installing Dask (:pr:`7308`) `Genevieve Buckley`_
- Support for Python 3.9 (:pr:`7289`) `Guido Imperiale`_
- Add broadcast-based merge implementation (:pr:`7143`) `Richard (Rick) Zamora`_
- Add ``split_every`` to ``graph_manipulation`` (:pr:`7282`) `Guido Imperiale`_
- Typo in optimize docs (:pr:`7306`) `Julius Busecke`_
- ``dask.graph_manipulation`` support for ``xarray.Dataset`` (:pr:`7276`) `Guido Imperiale`_
- Add plot width and height support for Bokeh 2.3.0 (:pr:`7297`) `James Bourbeau`_
- Add NumPy functions ``tri``, ``triu_indices``, ``triu_indices_from``, ``tril_indices``, ``tril_indices_from`` (:pr:`6997`) `Illviljan`_
- Remove "cleanup" task in DataFrame on-disk shuffle (:pr:`7260`) `Sinclair Target`_
- Use development version of ``distributed`` in CI (:pr:`7279`) `James Bourbeau`_
- Moving high level graph pack/unpack Dask  (:pr:`7179`) `Mads R. B. Kristensen`_
- Improve performance of ``merge_percentiles`` (:pr:`7172`) `Ashwin Srinath`_
- DOC: add ``dask-sql`` and ``fugue`` (:pr:`7129`) `Ray Bell`_
- Example for working with categoricals and parquet (:pr:`7085`) `McToel`_
- Adds tree reduction to ``bincount`` (:pr:`7183`) `Thomas J. Fan`_
- Improve documentation of ``name`` in ``from_array`` (:pr:`7264`) `Bruce Merry`_
- Fix ``cumsum`` for empty partitions (:pr:`7230`) `Julia Signell`_
- Add ``map_blocks`` example to dask array creation docs (:pr:`7221`) `Julia Signell`_
- Fix performance issue in ``dask.graph_manipulation.wait_on()`` (:pr:`7258`) `Guido Imperiale`_
- Replace coveralls with codecov.io (:pr:`7246`) `Guido Imperiale`_
- Pin to a particular ``black`` rev in pre-commit (:pr:`7256`) `Julia Signell`_
- Minor typo in documentation: ``array-chunks.rst`` (:pr:`7254`) `Magnus Nord`_
- Fix bugs in ``Blockwise`` and ``ShuffleLayer`` (:pr:`7213`) `Richard (Rick) Zamora`_
- Fix parquet filtering bug for ``"pyarrow-dataset"`` with pyarrow-3.0.0 (:pr:`7200`) `Richard (Rick) Zamora`_
- ``graph_manipulation`` without NumPy (:pr:`7243`) `Guido Imperiale`_
- Support for NEP-35 (:pr:`6738`) `Peter Andreas Entschev`_
- Avoid running unit tests during doctest CI build (:pr:`7240`) `James Bourbeau`_
- Run doctests on CI (:pr:`7238`) `Julia Signell`_
- Cleanup code quality on set arithmetics (:pr:`7196`) `Guido Imperiale`_
- Add ``dask.array.delete`` (:pr:`7125`) `Julia Signell`_
- Unpin graphviz now that new conda-forge recipe is built (:pr:`7235`) `Julia Signell`_
- Don't use NumPy 1.20 from conda-forge on Mac (:pr:`7211`) `Guido Imperiale`_
- ``map_overlap``: Don't rechunk axes without overlap (:pr:`7233`) `Deepak Cherian`_
- Pin graphviz to avoid issue with latest conda-forge build (:pr:`7232`) `Julia Signell`_
- Use ``html_css_files`` in docs for custom CSS (:pr:`7220`) `James Bourbeau`_
- Graph manipulation: ``clone``, ``bind``, ``checkpoint``, ``wait_on`` (:pr:`7109`) `Guido Imperiale`_
- Fix handling of filter expressions in parquet ``pyarrow-dataset`` engine (:pr:`7186`) `Joris Van den Bossche`_
- Extend ``__setitem__`` to more closely match numpy (:pr:`7033`) `David Hassell`_
- Clean up Python 2 syntax (:pr:`7195`) `Guido Imperiale`_
- Fix regression in ``Delayed._length`` (:pr:`7194`) `Guido Imperiale`_
- ``__dask_layers__()`` tests and tweaks (:pr:`7177`) `Guido Imperiale`_
- Properly convert ``HighLevelGraph`` in multiprocessing scheduler (:pr:`7191`) `Jim Crist-Harif`_
- Don't fail fast in CI (:pr:`7188`) `James Bourbeau`_


.. _v2021.02.0:

2021.02.0
---------

Released on February 5, 2021

- Add ``percentile`` support for NEP-35 (:pr:`7162`) `Peter Andreas Entschev`_
- Added support for ``Float64`` in column assignment (:pr:`7173`) `Nils Braun`_
- Coarsen rechunking error (:pr:`7127`) `Davis Bennett`_
- Fix upstream CI tests (:pr:`6896`) `Julia Signell`_
- Revise ``HighLevelGraph`` Mapping API (:pr:`7160`) `Guido Imperiale`_
- Update low-level graph spec to use any hashable for keys (:pr:`7163`) `James Bourbeau`_
- Generically rebuild a collection with different keys (:pr:`7142`) `Guido Imperiale`_
- Make easier to link issues in PRs (:pr:`7130`) `Ray Bell`_
- Add ``dask.array.append`` (:pr:`7146`) `D-Stacks`_
- Allow ``dask.array.ravel`` to accept ``array_like`` argument (:pr:`7138`) `D-Stacks`_
- Fixes link in array design doc (:pr:`7152`) `Thomas J. Fan`_
- Fix example of using ``blockwise`` for an outer product (:pr:`7119`) `Bruce Merry`_
- Deprecate ``HighlevelGraph.dicts`` in favor of ``.layers`` (:pr:`7145`) `Amit Kumar`_
- Align ``FastParquetEngine`` with pyarrow engines (:pr:`7091`) `Richard (Rick) Zamora`_
- Merge annotations (:pr:`7102`) `Ian Rose`_
- Simplify contents of parts list in ``read_parquet`` (:pr:`7066`) `Richard (Rick) Zamora`_
- ``check_meta(``): use ``__class__`` when checking DataFrame types (:pr:`7099`) `Mads R. B. Kristensen`_
- Cache several properties (:pr:`7104`) `Illviljan`_
- Fix parquet ``getitem`` optimization (:pr:`7106`) `Richard (Rick) Zamora`_
- Add cytoolz back to CI environment (:pr:`7103`) `James Bourbeau`_


.. _v2021.01.1:

2021.01.1
---------

Released on January 22, 2021

- Partially fix ``cumprod`` (:pr:`7089`) `Julia Signell`_
- Test pandas 1.1.x / 1.2.0 releases and pandas nightly (:pr:`6996`) `Joris Van den Bossche`_
- Use assign to avoid ``SettingWithCopyWarning`` (:pr:`7092`) `Julia Signell`_
- ``'mode'`` argument passed to ``bokeh.output_file()`` (:pr:`7034`) (:pr:`7075`) `patquem`_
- Skip empty partitions when doing ``groupby.value_counts`` (:pr:`7073`) `Julia Signell`_
- Add error messages to ``assert_eq()`` (:pr:`7083`) `James Lamb`_
- Make cached properties read-only (:pr:`7077`) `Illviljan`_


.. _v2021.01.0:

2021.01.0
---------

Released on January 15, 2021

- ``map_partitions`` with review comments (:pr:`6776`) `Kumar Bharath Prabhu`_
- Make sure that ``population`` is a real list (:pr:`7027`) `Julia Signell`_
- Propagate ``storage_options`` in ``read_csv`` (:pr:`7074`) `Richard (Rick) Zamora`_
- Remove all ``BlockwiseIO`` code (:pr:`7067`) `Richard (Rick) Zamora`_
- Fix CI (:pr:`7069`) `James Bourbeau`_
- Add option to control rechunking in ``reshape`` (:pr:`6753`) `Tom Augspurger`_
- Fix ``linalg.lstsq`` for complex inputs (:pr:`7056`) `Johnnie Gray`_
- Add ``compression='infer'`` default to ``read_csv`` (:pr:`6960`) `Richard (Rick) Zamora`_
- Revert parameter changes in ``svd_compressed`` #7003 (:pr:`7004`) `Eric Czech`_
- Skip failing s3 test (:pr:`7064`) `Martin Durant`_
- Revert ``BlockwiseIO`` (:pr:`7048`) `Richard (Rick) Zamora`_
- Add some cross-references to ``DataFrame.to_bag()`` and ``Series.to_bag()`` (:pr:`7049`) `Rob Malouf`_
- Rewrite ``matmul`` as ``blockwise`` without contraction/concatenate (:pr:`7000`) `Rafal Wojdyla`_
- Use ``functools.cached_property`` in ``da.shape`` (:pr:`7023`) `Illviljan`_
- Use meta value in series ``non_empty`` (:pr:`6976`) `Julia Signell`_
- Revert "Temporarly pin sphinx version to 3.3.1 (:pr:`7002`)" (:pr:`7014`) `Rafal Wojdyla`_
- Revert ``python-graphviz`` pinning (:pr:`7037`) `Julia Signell`_
- Accidentally committed print statement (:pr:`7038`) `Julia Signell`_
- Pass ``dropna`` and ``observed`` in ``agg`` (:pr:`6992`) `Julia Signell`_
- Add index to ``meta`` after ``.str.split`` with expand (:pr:`7026`) `Ruben van de Geer`_
- CI: test pyarrow 2.0 and nightly (:pr:`7030`) `Joris Van den Bossche`_
- Temporarily pin ``python-graphviz`` in CI (:pr:`7031`) `James Bourbeau`_
- Underline section in ``numpydoc`` (:pr:`7013`) `Matthias Bussonnier`_
- Keep normal optimizations when adding custom optimizations (:pr:`7016`) `Matthew Rocklin`_
- Temporarily pin sphinx version to 3.3.1 (:pr:`7002`) `Rafal Wojdyla`_
- DOC: Misc formatting (:pr:`6998`) `Matthias Bussonnier`_
- Add ``inline_array`` option to ``from_array`` (:pr:`6773`) `Tom Augspurger`_
- Revert "Initial pass at blockwise array creation routines (:pr:`6931)" (:pr:`6995`) `James Bourbeau`_
- Set ``npartitions`` in ``set_index`` (:pr:`6978`) `Julia Signell`_
- Upstream ``config`` serialization and inheritance (:pr:`6987`) `Jacob Tomlinson`_
- Bump the minimum time in ``test_minimum_time`` (:pr:`6988`) `Martin Durant`_
- Fix pandas ``dtype`` inference for ``read_parquet`` (:pr:`6985`) `Richard (Rick) Zamora`_
- Avoid data loss in ``set_index`` with ``sorted=True`` (:pr:`6980`) `Richard (Rick) Zamora`_
- Bugfix in ``read_parquet`` for handling un-named indices with ``index=False`` (:pr:`6969`) `Richard (Rick) Zamora`_
- Use ``__class__`` when comparing meta data (:pr:`6981`) `Mads R. B. Kristensen`_
- Comparing string versions won't always work (:pr:`6979`) `Rafal Wojdyla`_
- Fix :pr:`6925` (:pr:`6982`) `sdementen`_
- Initial pass at blockwise array creation routines (:pr:`6931`) `Ian Rose`_
- Simplify ``has_parallel_type()`` (:pr:`6927`) `Mads R. B. Kristensen`_
- Handle annotation unpacking in ``BlockwiseIO`` (:pr:`6934`) `Simon Perkins`_
- Avoid deprecated ``yield_fixture`` in ``test_sql.py`` (:pr:`6968`) `Richard (Rick) Zamora`_
- Remove bad graph logic in ``BlockwiseIO`` (:pr:`6933`) `Richard (Rick) Zamora`_
- Get config item if variable is ``None`` (:pr:`6862`) `Jacob Tomlinson`_
- Update ``from_pandas`` docstring (:pr:`6957`) `Richard (Rick) Zamora`_
- Prevent ``fuse_roots`` from clobbering annotations (:pr:`6955`) `Simon Perkins`_


.. _v2020.12.0:

2020.12.0
---------

Released on December 10, 2020

Highlights
^^^^^^^^^^

- Switched to `CalVer <https://calver.org/>`_ for versioning scheme.
- Introduced new APIs for ``HighLevelGraph`` to enable sending high-level representations of
  task graphs to the distributed scheduler.
- Introduced new ``HighLevelGraph`` layer objects including ``BasicLayer``, ``Blockwise``,
  ``BlockwiseIO``, ``ShuffleLayer``, and more.
- Added support for applying custom ``Layer``-level annotations like ``priority``, ``retries``,
  etc. with the ``dask.annotations`` context manager.
- Updated minimum supported version of pandas to 0.25.0 and NumPy to 1.15.1.
- Support for the ``pyarrow.dataset`` API to ``read_parquet``.
- Several fixes to Dask Array's SVD.

All changes
^^^^^^^^^^^

- Make ``observed`` kwarg optional (:pr:`6952`) `Julia Signell`_
- Min supported pandas 0.25.0 numpy 1.15.1 (:pr:`6895`) `Julia Signell`_
- Make order of categoricals unambiguous (:pr:`6949`) `Julia Signell`_
- Improve "pyarrow-dataset" statistics performance for ``read_parquet`` (:pr:`6918`) `Richard (Rick) Zamora`_
- Add ``observed`` keyword to ``groupby`` (:pr:`6854`) `Julia Signell`_
- Make sure ``include_path_column`` works when there are multiple partitions per file (:pr:`6911`) `Julia Signell`_
- Fix: ``array.overlap`` and ``array.map_overlap`` block sizes are incorrect when depth is an unsigned bit type (:pr:`6909`) `GFleishman`_
- Fix syntax error in HLG docs example (:pr:`6946`) `Mark`_
- Return a ``Bag`` from ``sample`` (:pr:`6941`) `Shang Wang`_
- Add ``ravel_multi_index`` (:pr:`6939`) `Illviljan`_
- Enable parquet metadata collection in parallel (:pr:`6921`) `Richard (Rick) Zamora`_
- Avoid using ``_file`` in ``progressbar`` if it is ``None`` (:pr:`6938`) `Mark Harfouche`_
- Add Zarr to upstream CI build (:pr:`6932`) `James Bourbeau`_
- Introduce ``BlockwiseIO`` layer (:pr:`6878`) `Richard (Rick) Zamora`_
- Transmit ``Layer`` Annotations to Scheduler (:pr:`6889`) `Simon Perkins`_
- Update opportunistic caching page to remove experimental warning (:pr:`6926`) `Timost`_
- Allow ``pyarrow >2.0.0`` (:pr:`6772`) `Richard (Rick) Zamora`_
- Support ``pyarrow.dataset`` API for ``read_parquet`` (:pr:`6534`) `Richard (Rick) Zamora`_
- Add more informative error message to ``da.coarsen`` when coarsening factors do not divide shape (:pr:`6908`) `Davis Bennett`_
- Only run the cron CI on ``dask/dask`` not forks (:pr:`6905`) `Jacob Tomlinson`_
- Add ``annotations`` to ``ShuffleLayers`` (:pr:`6913`) `Matthew Rocklin`_
- Temporarily xfail ``test_from_s3`` (:pr:`6915`) `James Bourbeau`_
- Added dataframe ``skew`` method (:pr:`6881`) `Jan Borchmann`_
- Fix ``dtype`` in array ``meta`` (:pr:`6893`) `Julia Signell`_
- Missing ``name`` arg in ``helm install ...`` (:pr:`6903`) `Ruben van de Geer`_
- Fix: exception when reading an item with filters (:pr:`6901`) `Martin Durant`_
- Add support for ``cupyx`` sparse to ``dask.array.dot`` (:pr:`6846`) `Akira Naruse`_
- Pin array mindeps up a bit to get the tests to pass [test-mindeps] (:pr:`6894`) `Julia Signell`_
- Update/remove pandas and numpy in mindeps (:pr:`6888`) `Julia Signell`_
- Fix ``ArrowEngine`` bug in use of ``clear_known_categories`` (:pr:`6887`) `Richard (Rick) Zamora`_
- Fix documentation about task scheduler (:pr:`6879`) `Zhengnan Zhao`_
- Add human relative time formatting utility (:pr:`6883`) `Jacob Tomlinson`_
- Possible fix for 6864 ``set_index`` issue (:pr:`6866`) `Richard (Rick) Zamora`_
- ``BasicLayer``: remove dependency arguments (:pr:`6859`) `Mads R. B. Kristensen`_
- Serialization of ``Blockwise`` (:pr:`6848`) `Mads R. B. Kristensen`_
- Address ``columns=[]`` bug (:pr:`6871`) `Richard (Rick) Zamora`_
- Avoid duplicate parquet schema communication (:pr:`6841`) `Richard (Rick) Zamora`_
- Add ``create_metadata_file`` utility for existing parquet datasets (:pr:`6851`) `Richard (Rick) Zamora`_
- Improve ordering for workloads with a common terminus (:pr:`6779`) `Tom Augspurger`_
- Stringify utilities (:pr:`6852`) `Mads R. B. Kristensen`_
- Add keyword ``overwrite=True`` to ``to_parquet`` to remove dangling files when overwriting a pyarrow ``Dataset``. (:pr:`6825`) `Greg Hayes`_
- Removed ``map_tasks()`` and ``map_basic_layers()`` (:pr:`6853`) `Mads R. B. Kristensen`_
- Introduce QR iteration to ``svd_compressed`` (:pr:`6813`) `RogerMoens`_
- ``__dask_distributed_pack__()`` now takes a ``client`` argument (:pr:`6850`) `Mads R. B. Kristensen`_
- Use ``map_partitions`` instead of ``delayed`` in ``set_index`` (:pr:`6837`) `Mads R. B. Kristensen`_
- Add doc hit for ``as_completed().update(futures)`` (:pr:`6817`) `manuels`_
- Bump GHA ``setup-miniconda`` version (:pr:`6847`) `Jacob Tomlinson`_
- Remove nans when setting sorted index (:pr:`6829`) `Rockwell Weiner`_
- Fix transpose of u in SVD (:pr:`6799`) `RogerMoens`_
- Migrate to GitHub Actions (:pr:`6794`) `Jacob Tomlinson`_
- Fix sphinx ``currentmodule`` usage (:pr:`6839`) `James Bourbeau`_
- Fix minimum dependencies CI builds (:pr:`6838`) `James Bourbeau`_
- Avoid graph materialization during ``Blockwise`` culling (:pr:`6815`) `Richard (Rick) Zamora`_
- Fixed typo (:pr:`6834`) `Devanshu Desai`_
- Use ``HighLevelGraph.merge`` in ``collections_to_dsk`` (:pr:`6836`) `Mads R. B. Kristensen`_
- Respect ``dtype`` in svd ``compression_matrix`` #2849 (:pr:`6802`) `RogerMoens`_
- Add blocksize to task name (:pr:`6818`) `Julia Signell`_
- Check for all-NaN partitions (:pr:`6821`) `Rockwell Weiner`_
- Change "institutional" SQL doc section to point to main SQL doc (:pr:`6823`) `Martin Durant`_
- Fix: ``DataFrame.join`` doesn't accept Series as other (:pr:`6809`) `David Katz`_
- Remove ``to_delayed`` operations from ``to_parquet`` (:pr:`6801`) `Richard (Rick) Zamora`_
- Layer annotation docstrings improvements (:pr:`6806`) `Simon Perkins`_
- Avro reader (:pr:`6780`) `Martin Durant`_
- Rechunk array if smallest chunk size is smaller than depth (:pr:`6708`) `Julia Signell`_
- Add Layer Annotations (:pr:`6767`) `Simon Perkins`_
- Add "view code" links to documentation (:pr:`6793`) `manuels`_
- Add optional IO-subgraph to ``Blockwise`` Layers (:pr:`6715`) `Richard (Rick) Zamora`_
- Add high level graph pack/unpack for distributed (:pr:`6786`) `Mads R. B. Kristensen`_
- Add missing methods of the Dataframe API (:pr:`6789`) `Stephannie Jimenez Gacha`_
- Add doc on managing environments (:pr:`6778`) `Martin Durant`_
- HLG: ``get_all_external_keys()`` (:pr:`6774`) `Mads R. B. Kristensen`_
- Avoid rechunking in reshape with ``chunksize=1`` (:pr:`6748`) `Tom Augspurger`_
- Try to make categoricals work on join (:pr:`6205`) `Julia Signell`_
- Fix some minor typos and trailing whitespaces in ``array-slice.rst`` (:pr:`6771`) `Magnus Nord`_
- Bugfix for parquet metadata writes of empty dataframe partitions (pyarrow)  (:pr:`6741`) `Callum Noble`_
- Document ``meta`` kwarg in ``map_blocks`` and ``map_overlap``. (:pr:`6763`) `Peter Andreas Entschev`_
- Begin experimenting with parallel prefix scan for ``cumsum`` and ``cumprod`` (:pr:`6675`) `Erik Welch`_
- Clarify differences in boolean indexing between dask and numpy arrays (:pr:`6764`) `Illviljan`_
- Efficient serialization of shuffle layers (:pr:`6760`) `James Bourbeau`_
- Config array optimize to skip fusion and return a HLG (:pr:`6751`) `Mads R. B. Kristensen`_
- Temporarily use ``pyarrow<2`` in CI (:pr:`6759`) `James Bourbeau`_
- Fix meta for ``min``/``max`` reductions (:pr:`6736`) `Peter Andreas Entschev`_
- Add 2D possibility to ``da.linalg.lstsq`` - mirroring numpy (:pr:`6749`) `Pascal Bourgault`_
- CI: Fixed bug causing flaky test failure in pivot (:pr:`6752`) `Tom Augspurger`_
- Serialization of layers (:pr:`6693`) `Mads R. B. Kristensen`_
- Add ``attrs`` property to Series/Dataframe (:pr:`6742`) `Illviljan`_
- Removed Mutable Default Argument (:pr:`6747`) `Mads R. B. Kristensen`_
- Adjust parquet ``ArrowEngine`` to allow more easy subclass for writing (:pr:`6505`) `Joris Van den Bossche`_
- Add ``ShuffleStage`` HLG Layer (:pr:`6650`) `Richard (Rick) Zamora`_
- Handle literal in ``meta_from_array`` (:pr:`6731`) `Peter Andreas Entschev`_
- Do balanced rechunking even if chunks are the same (:pr:`6735`) `Chris Roat`_
- Fix docstring ``DataFrame.set_index`` (:pr:`6739`) `Gil Forsyth`_
- Ensure ``HighLevelGraph`` layers always contain ``Layer`` instances (:pr:`6716`) `James Bourbeau`_
- Map on ``HighLevelGraph`` Layers (:pr:`6689`) `Mads R. B. Kristensen`_
- Update overlap ``*_like`` function calls and CuPy tests (:pr:`6728`) `Peter Andreas Entschev`_
- Fixes for ``svd`` with ``__array_function__`` (:pr:`6727`) `Peter Andreas Entschev`_
- Added doctest extension for documentation (:pr:`6397`) `Jim Circadian`_
- Minor fix to #5628 using @pentschev's suggestion (:pr:`6724`) `John A Kirkham`_
- Change type of Dask array when meta type changes (:pr:`5628`) `Matthew Rocklin`_
- Add ``az`` (:pr:`6719`) `Ray Bell`_
- HLG: ``get_dependencies()`` of single keys (:pr:`6699`) `Mads R. B. Kristensen`_
- Revert "Revert "Use HighLevelGraph layers everywhere in collections (:pr:`6510`)" (:pr:`6697`)" (:pr:`6707`) `Tom Augspurger`_
- Allow ``*_like`` array creation functions to respect input array type (:pr:`6680`) `Genevieve Buckley`_
- Update ``dask-sphinx-theme`` version (:pr:`6700`) `Gil Forsyth`_


.. _v2.30.0 / 2020-10-06:

2.30.0 / 2020-10-06
-------------------

Array
^^^^^

- Allow ``rechunk`` to evenly split into N chunks (:pr:`6420`) `Scott Sievert`_


.. _v2.29.0 / 2020-10-02:

2.29.0 / 2020-10-02
-------------------

Array
^^^^^

- ``_repr_html_``: color sides darker instead of drawing all the lines (:pr:`6683`) `Julia Signell`_
- Removes warning from ``nanstd`` and ``nanvar`` (:pr:`6667`) `Thomas J. Fan`_
- Get shape of output from original array - ``map_overlap`` (:pr:`6682`) `Julia Signell`_
- Replace ``np.searchsorted`` with ``bisect`` in indexing (:pr:`6669`) `Joachim B Haga`_

Bag
^^^

- Make sure subprocesses have a consistent hash for bag ``groupby`` (:pr:`6660`) `Itamar Turner-Trauring`_

Core
^^^^

- Revert "Use ``HighLevelGraph`` layers everywhere in collections (:pr:`6510`)" (:pr:`6697`) `Tom Augspurger`_
- Use ``pandas.testing`` (:pr:`6687`) `John A Kirkham`_
- Improve 128-bit floating-point skip in tests (:pr:`6676`) `Elliott Sales de Andrade`_

DataFrame
^^^^^^^^^

- Allow setting dataframe items using a bool dataframe (:pr:`6608`) `Julia Signell`_

Documentation
^^^^^^^^^^^^^

- Fix typo (:pr:`6692`) `garanews`_
- Fix a few typos (:pr:`6678`) `Pav A`_


.. _v2.28.0 / 2020-09-25:

2.28.0 / 2020-09-25
-------------------

Array
^^^^^

- Partially reverted changes to ``Array`` indexing that produces large changes.
  This restores the behavior from Dask 2.25.0 and earlier, with a warning
  when large chunks are produced. A configuration option is provided
  to avoid creating the large chunks, see :ref:`array.slicing.efficiency`.
  (:pr:`6665`) `Tom Augspurger`_
- Add ``meta`` to ``to_dask_array`` (:pr:`6651`) `Kyle Nicholson`_
- Fix :pr:`6631` and :pr:`6611` (:pr:`6632`) `Rafal Wojdyla`_
- Infer object in array reductions (:pr:`6629`) `Daniel Saxton`_
- Adding ``v_based`` flag for ``svd_flip`` (:pr:`6658`) `Eric Czech`_
- Fix flakey array ``mean`` (:pr:`6656`) `Sam Grayson`_

Core
^^^^

- Removed ``dsk`` equality check from ``SubgraphCallable.__eq__`` (:pr:`6666`) `Mads R. B. Kristensen`_
- Use ``HighLevelGraph`` layers everywhere in collections (:pr:`6510`) `Mads R. B. Kristensen`_
- Adds hash dunder method to ``SubgraphCallable`` for caching purposes (:pr:`6424`) `Andrew Fulton`_
- Stop writing commented out config files by default (:pr:`6647`) `Matthew Rocklin`_

DataFrame
^^^^^^^^^

- Add support for collect list aggregation via ``agg`` API (:pr:`6655`) `Madhur Tandon`_
- Slightly better error message (:pr:`6657`) `Julia Signell`_


.. _v2.27.0 / 2020-09-18:

2.27.0 / 2020-09-18
-------------------

Array
^^^^^

- Preserve ``dtype`` in ``svd`` (:pr:`6643`) `Eric Czech`_

Core
^^^^

- ``store()``: create a single HLG layer (:pr:`6601`) `Mads R. B. Kristensen`_
- Add pre-commit CI build (:pr:`6645`) `James Bourbeau`_
- Update ``.pre-commit-config`` to latest black. (:pr:`6641`) `Julia Signell`_
- Update super usage to remove Python 2 compatibility (:pr:`6630`) `Poruri Sai Rahul`_
- Remove u string prefixes (:pr:`6633`) `Poruri Sai Rahul`_

DataFrame
^^^^^^^^^

- Improve error message for ``to_sql`` (:pr:`6638`) `Julia Signell`_
- Use empty list as categories (:pr:`6626`) `Julia Signell`_

Documentation
^^^^^^^^^^^^^

- Add ``autofunction`` to array api docs for more ufuncs (:pr:`6644`) `James Bourbeau`_
- Add a number of missing ufuncs to ``dask.array`` docs (:pr:`6642`) `Ralf Gommers`_
- Add ``HelmCluster`` docs (:pr:`6290`) `Jacob Tomlinson`_


.. _v2.26.0 / 2020-09-11:

2.26.0 / 2020-09-11
-------------------

Array
^^^^^

- Backend-aware dtype inference for single-chunk svd (:pr:`6623`) `Eric Czech`_
- Make ``array.reduction`` docstring match for dtype (:pr:`6624`) `Martin Durant`_
- Set lower bound on compression level for ``svd_compressed`` using rows and cols (:pr:`6622`) `Eric Czech`_
- Improve SVD consistency and small array handling (:pr:`6616`) `Eric Czech`_
- Add ``svd_flip`` #6599 (:pr:`6613`) `Eric Czech`_
- Handle sequences containing dask Arrays (:pr:`6595`) `Gabe Joseph`_
- Avoid large chunks from ``getitem`` with lists (:pr:`6514`) `Tom Augspurger`_
- Eagerly slice numpy arrays in ``from_array`` (:pr:`6605`) `Deepak Cherian`_
- Restore ability to pickle dask arrays (:pr:`6594`) `Noah D. Brenowitz`_
- Add SVD support for short-and-fat arrays (:pr:`6591`) `Eric Czech`_
- Add simple chunk type registry and defer as appropriate to upcast types (:pr:`6393`) `Jon Thielen`_
- Align coarsen chunks by default (:pr:`6580`) `Deepak Cherian`_
- Fixup reshape on unknown dimensions and other testing fixes (:pr:`6578`) `Ryan Williams`_

Core
^^^^

- Add validation and fixes for ``HighLevelGraph`` dependencies (:pr:`6588`) `Mads R. B. Kristensen`_
- Fix linting issue (:pr:`6598`) `Tom Augspurger`_
- Skip ``bokeh`` version 2.0.0 (:pr:`6572`) `John A Kirkham`_

DataFrame
^^^^^^^^^

- Added bytes/row calculation when using meta (:pr:`6585`) `McToel`_
- Handle ``min_count`` in ``Series.sum`` / ``prod`` (:pr:`6618`) `Daniel Saxton`_
- Update ``DataFrame.set_index`` docstring (:pr:`6549`) `Timost`_
- Always compute 0 and 1 quantiles during quantile calculations (:pr:`6564`) `Erik Welch`_
- Fix wrong path when reading empty csv file (:pr:`6573`) `Abdulelah Bin Mahfoodh`_

Documentation
^^^^^^^^^^^^^

- Doc: Troubleshooting dashboard 404 (:pr:`6215`) `Kilian Lieret`_
- Fixup ``extraConfig`` example (:pr:`6625`) `Tom Augspurger`_
- Update supported Python versions (:pr:`6609`) `Julia Signell`_
- Document dask/daskhub helm chart (:pr:`6560`) `Tom Augspurger`_


.. _v2.25.0 / 2020-08-28:

2.25.0 / 2020-08-28
-------------------

Core
^^^^

- Compare key hashes in ``subs()`` (:pr:`6559`) `Mads R. B. Kristensen`_
- Rerun with latest ``black`` release (:pr:`6568`) `James Bourbeau`_
- License update (:pr:`6554`) `Tom Augspurger`_

DataFrame
^^^^^^^^^

- Add gs ``read_parquet`` example (:pr:`6548`) `Ray Bell`_

Documentation
^^^^^^^^^^^^^

- Remove version from documentation page names (:pr:`6558`) `James Bourbeau`_
- Update ``kubernetes-helm.rst`` (:pr:`6523`) `David Sheldon`_
- Stop 2020 survey (:pr:`6547`) `Tom Augspurger`_


.. _v2.24.0 / 2020-08-22:

2.24.0 / 2020-08-22
-------------------

Array
^^^^^

-   Fix setting random seed in tests. (:pr:`6518`) `Elliott Sales de Andrade`_
-   Support meta in apply gufunc (:pr:`6521`) `joshreback`_
-   Replace `cupy.sparse` with `cupyx.scipy.sparse` (:pr:`6530`) `John A Kirkham`_

Dataframe
^^^^^^^^^

-   Bump up tolerance for rolling tests (:pr:`6502`) `Julia Signell`_
-   Implement DatFrame.__len__ (:pr:`6515`) `Tom Augspurger`_
-   Infer arrow schema in to_parquet  (for ArrowEngine`) (:pr:`6490`) `Richard (Rick) Zamora`_
-   Fix parquet test when no pyarrow (:pr:`6524`) `Martin Durant`_
-   Remove problematic ``filter`` arguments in ArrowEngine (:pr:`6527`) `Richard (Rick) Zamora`_
-   Avoid schema validation by default in ArrowEngine (:pr:`6536`) `Richard (Rick) Zamora`_

Core
^^^^

-   Use unpack_collections in make_blockwise_graph (:pr:`6517`) `Thomas J. Fan`_
-   Move key_split() from optimization.py to utils.py (:pr:`6529`) `Mads R. B. Kristensen`_
-   Make tests run on moto server (:pr:`6528`) `Martin Durant`_


.. _v2.23.0 / 2020-08-14:

2.23.0 / 2020-08-14
-------------------

Array
^^^^^

- Reduce ``np.zeros``, ``ones``, and ``full`` array size with broadcasting (:pr:`6491`) `Matthias Bussonnier`_
- Add missing ``meta=`` for ``trim`` in ``map_overlap`` (:pr:`6494`) `Peter Andreas Entschev`_

Bag
^^^

- Bag repartition partition size (:pr:`6371`) `joshreback`_

Core
^^^^

- ``Scalar.__dask_layers__()`` to return ``self._name`` instead of ``self.key`` (:pr:`6507`) `Mads R. B. Kristensen`_
- Update dependencies correctly in ``fuse_root`` optimization (:pr:`6508`) `Mads R. B. Kristensen`_


DataFrame
^^^^^^^^^

- Adds ``items`` to dataframe (:pr:`6503`) `Thomas J. Fan`_
- Include compression in ``write_table`` call (:pr:`6499`) `Julia Signell`_
- Fixed warning in ``nonempty_series`` (:pr:`6485`) `Tom Augspurger`_
- Intelligently determine partitions based on type of first arg (:pr:`6479`) `Matthew Rocklin`_
- Fix pyarrow ``mkdirs`` (:pr:`6475`) `Julia Signell`_
- Fix duplicate parquet output in ``to_parquet`` (:pr:`6451`) `michaelnarodovitch`_

Documentation
^^^^^^^^^^^^^

- Fix documentation ``da.histogram`` (:pr:`6439`) `Roberto Panai`_
- Add ``agg`` ``nunique`` example (:pr:`6404`) `Ray Bell`_
- Fixed a few typos in the SQL docs (:pr:`6489`) `Mike McCarty`_
- Docs for SQLing (:pr:`6453`) `Martin Durant`_


.. _v2.22.0 / 2020-07-31:

2.22.0 / 2020-07-31
-------------------

Array
^^^^^

- Compatibility for NumPy dtype deprecation (:pr:`6430`) `Tom Augspurger`_

Core
^^^^

- Implement ``sizeof`` for some ``bytes``-like objects (:pr:`6457`) `John A Kirkham`_
- HTTP error for new ``fsspec`` (:pr:`6446`) `Martin Durant`_
- When ``RecursionError`` is raised, return uuid from ``tokenize`` function (:pr:`6437`) `Julia Signell`_
- Install deps of upstream-dev packages (:pr:`6431`) `Tom Augspurger`_
- Use updated link in ``setup.cfg`` (:pr:`6426`) `Zhengnan Zhao`_

DataFrame
^^^^^^^^^

- Add single quotes around column names if strings (:pr:`6471`) `Gil Forsyth`_
- Refactor ``ArrowEngine`` for better ``read_parquet`` performance (:pr:`6346`) `Richard (Rick) Zamora`_
- Add ``tolist`` dispatch (:pr:`6444`) `GALI PREM SAGAR`_
- Compatibility with pandas 1.1.0rc0 (:pr:`6429`) `Tom Augspurger`_
- Multi value pivot table (:pr:`6428`) `joshreback`_
- Duplicate argument definitions in ``to_csv`` docstring (:pr:`6411`) `Jun Han (Johnson) Ooi`_

Documentation
^^^^^^^^^^^^^

- Add utility to docs to convert YAML config to env vars and back (:pr:`6472`) `Jacob Tomlinson`_
- Fix parameter server rendering (:pr:`6466`) `Scott Sievert`_
- Fixes broken links (:pr:`6403`) `Jim Circadian`_
- Complete parameter server implementation in docs (:pr:`6449`) `Scott Sievert`_
- Fix typo (:pr:`6436`) `Jack Xiaosong Xu`_


.. _v2.21.0 / 2020-07-17:

2.21.0 / 2020-07-17
-------------------

Array
^^^^^

- Correct error message in ``array.routines.gradient()`` (:pr:`6417`) `johnomotani`_
- Fix blockwise concatenate for array with some ``dimension=1`` (:pr:`6342`) `Matthias Bussonnier`_

Bag
^^^

- Fix ``bag.take`` example (:pr:`6418`) `Roberto Panai`_

Core
^^^^

- Groups values in optimization pass should only be graph and keys -- not an optimization + keys (:pr:`6409`) `Benjamin Zaitlen`_
- Call custom optimizations once, with ``kwargs`` provided (:pr:`6382`) `Clark Zinzow`_
- Include ``pickle5`` for testing on Python 3.7 (:pr:`6379`) `John A Kirkham`_

DataFrame
^^^^^^^^^

- Correct typo in error message (:pr:`6422`) `Tom McTiernan`_
- Use ``pytest.warns`` to check for ``UserWarning`` (:pr:`6378`) `Richard (Rick) Zamora`_
- Parse ``bytes_per_chunk keyword`` from string (:pr:`6370`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^

- Numpydoc formatting (:pr:`6421`) `Matthias Bussonnier`_
- Unpin ``numpydoc`` following 1.1 release (:pr:`6407`) `Gil Forsyth`_
- Numpydoc formatting (:pr:`6402`) `Matthias Bussonnier`_
- Add instructions for using conda when installing code for development (:pr:`6399`) `Ray Bell`_
- Update ``visualize`` docstrings (:pr:`6383`) `Zhengnan Zhao`_


.. _v2.20.0 / 2020-07-02:

2.20.0 / 2020-07-02
-------------------

Array
^^^^^

- Register ``sizeof`` for numpy zero-strided arrays (:pr:`6343`) `Matthias Bussonnier`_
- Use ``concatenate_lookup`` in ``concatenate`` (:pr:`6339`) `John A Kirkham`_
- Fix rechunking of arrays with some zero-length dimensions (:pr:`6335`) `Matthias Bussonnier`_

DataFrame
^^^^^^^^^

- Dispatch ``iloc``` calls to ``getitem`` (:pr:`6355`) `Gil Forsyth`_
- Handle unnamed pandas ``RangeIndex`` in fastparquet engine (:pr:`6350`) `Richard (Rick) Zamora`_
- Preserve index when writing partitioned parquet datasets with pyarrow (:pr:`6282`) `Richard (Rick) Zamora`_
- Use ``ignore_index`` for pandas' ``group_split_dispatch`` (:pr:`6251`) `Richard (Rick) Zamora`_

Documentation
^^^^^^^^^^^^^

- Add doc describing argument (:pr:`6318`) `asmith26`_


.. _v2.19.0 / 2020-06-19:

2.19.0 / 2020-06-19
-------------------

Array
^^^^^

- Cast chunk sizes to python int ``dtype`` (:pr:`6326`) `Gil Forsyth`_
- Add ``shape=None`` to ``*_like()`` array creation functions (:pr:`6064`) `Anderson Banihirwe`_

Core
^^^^

- Update expected error msg for protocol difference in fsspec (:pr:`6331`) `Gil Forsyth`_
- Fix for floats < 1 in ``parse_bytes`` (:pr:`6311`) `Gil Forsyth`_
- Fix exception causes all over the codebase (:pr:`6308`) `Ram Rachum`_
- Fix duplicated tests (:pr:`6303`) `James Lamb`_
- Remove unused testing function (:pr:`6304`) `James Lamb`_

DataFrame
^^^^^^^^^

- Add high-level CSV Subgraph (:pr:`6262`) `Gil Forsyth`_
- Fix ``ValueError`` when merging an index-only 1-partition dataframe (:pr:`6309`) `Krishan Bhasin`_
- Make ``index.map`` clear divisions. (:pr:`6285`) `Julia Signell`_

Documentation
^^^^^^^^^^^^^

- Add link to 2020 survey (:pr:`6328`) `Tom Augspurger`_
- Update ``bag.rst`` (:pr:`6317`) `Ben Shaver`_


.. _v2.18.1 / 2020-06-09:

2.18.1 / 2020-06-09
-------------------

Array
^^^^^

- Don't try to set name on ``full`` (:pr:`6299`) `Julia Signell`_
- Histogram: support lazy values for range/bins (another way) (:pr:`6252`) `Gabe Joseph`_

Core
^^^^

- Fix exception causes in ``utils.py`` (:pr:`6302`) `Ram Rachum`_
- Improve performance of ``HighLevelGraph`` construction (:pr:`6293`) `Julia Signell`_

Documentation
^^^^^^^^^^^^^

- Now readthedocs builds unrelased features' docstrings (:pr:`6295`) `Antonio Ercole De Luca`_
- Add ``asyncssh`` intersphinx mappings (:pr:`6298`) `Jacob Tomlinson`_


.. _v2.18.0 / 2020-06-05:

2.18.0 / 2020-06-05
-------------------

Array
^^^^^

- Cast slicing index to dask array if same shape as original (:pr:`6273`) `Julia Signell`_
- Fix ``stack`` error message (:pr:`6268`) `Stephanie Gott`_
- ``full`` & ``full_like``: error on non-scalar ``fill_value`` (:pr:`6129`) `Huite`_
- Support for multiple arrays in ``map_overlap`` (:pr:`6165`) `Eric Czech`_
- Pad resample divisions so that edges are counted (:pr:`6255`) `Julia Signell`_

Bag
^^^

- Random sampling of k elements from a dask bag #4799 (:pr:`6239`) `Antonio Ercole De Luca`_

DataFrame
^^^^^^^^^

- Add ``dropna``, ``sort``, and ``ascending`` to ``sort_values`` (:pr:`5880`) `Julia Signell`_
- Generalize ``from_dask_array`` (:pr:`6263`) `GALI PREM SAGAR`_
- Add derived docstring for ``SeriesGroupby.nunique`` (:pr:`6284`) `Julia Signell`_
- Remove ``NotImplementedError`` in resample with rule  (:pr:`6274`) `Abdulelah Bin Mahfoodh`_
- Add ``dd.to_sql`` (:pr:`6038`) `Ryan Williams`_

Documentation
^^^^^^^^^^^^^

- Update remote data section (:pr:`6258`) `Ray Bell`_


.. _v2.17.2 / 2020-05-28:

2.17.2 / 2020-05-28
-------------------

Core
^^^^

- Re-add the ``complete`` extra (:pr:`6257`) `Jim Crist-Harif`_

DataFrame
^^^^^^^^^

- Raise error if ``resample`` isn't going to give right answer (:pr:`6244`) `Julia Signell`_


.. _v2.17.1 / 2020-05-28:

2.17.1 / 2020-05-28
-------------------

Array
^^^^^

- Empty array rechunk (:pr:`6233`) `Andrew Fulton`_

Core
^^^^

- Make ``pyyaml`` required (:pr:`6250`) `Jim Crist-Harif`_
- Fix install commands from ``ImportError`` (:pr:`6238`) `Gaurav Sheni`_
- Remove issue template (:pr:`6249`) `Jacob Tomlinson`_

DataFrame
^^^^^^^^^

- Pass ``ignore_index`` to ``dd_shuffle`` from ``DataFrame.shuffle`` (:pr:`6247`) `Richard (Rick) Zamora`_
- Cope with missing HDF keys (:pr:`6204`) `Martin Durant`_
- Generalize ``describe`` & ``quantile`` apis (:pr:`5137`) `GALI PREM SAGAR`_


.. _v2.17.0 / 2020-05-26:

2.17.0 / 2020-05-26
-------------------

Array
^^^^^

- Small improvements to ``da.pad`` (:pr:`6213`) `Mark Boer`_
- Return ``tuple`` if multiple outputs in ``dask.array.apply_gufunc``, add test to check for tuple (:pr:`6207`) `Kai Mühlbauer`_
- Support ``stack`` with unknown chunksizes (:pr:`6195`) `swapna`_

Bag
^^^

- Random Choice on Bags (:pr:`6208`) `Antonio Ercole De Luca`_

Core
^^^^

- Raise warning ``delayed.visualise()`` (:pr:`6216`) `Amol Umbarkar`_
- Ensure other pickle arguments work (:pr:`6229`) `John A Kirkham`_
- Overhaul ``fuse()`` config (:pr:`6198`) `Guido Imperiale`_
- Update ``dask.order.order`` to consider "next" nodes using both FIFO and LIFO (:pr:`5872`) `Erik Welch`_

DataFrame
^^^^^^^^^

- Use 0 as ``fill_value`` for more agg methods (:pr:`6245`) `Julia Signell`_
- Generalize ``rearrange_by_column_tasks`` and add ``DataFrame.shuffle`` (:pr:`6066`) `Richard (Rick) Zamora`_
- Xfail ``test_rolling_numba_engine`` for newer numba and older pandas (:pr:`6236`) `James Bourbeau`_
- Generalize ``fix_overlap`` (:pr:`6240`) `GALI PREM SAGAR`_
- Fix ``DataFrame.shape`` with no columns (:pr:`6237`) `noreentry`_
- Avoid shuffle when setting a presorted index with overlapping divisions (:pr:`6226`) `Krishan Bhasin`_
- Adjust the Parquet engine classes to allow more easily subclassing (:pr:`6211`) `Marius van Niekerk`_
- Fix ``dd.merge_asof`` with ``left_on='col'`` & ``right_index=True`` (:pr:`6192`) `noreentry`_
- Disable warning for ``concat`` (:pr:`6210`) `Tung Dang`_
- Move ``AUTO_BLOCKSIZE`` out of ``read_csv`` signature (:pr:`6214`) `Jim Crist-Harif`_
- ``.loc`` indexing with callable (:pr:`6185`) `Endre Mark Borza`_
- Avoid apply in ``_compute_sum_of_squares`` for groupby std agg (:pr:`6186`) `Richard (Rick) Zamora`_
- Minor correction to ``test_parquet`` (:pr:`6190`) `Brian Larsen`_
- Adhering to the passed pat for delimeter join and fix error message (:pr:`6194`) `GALI PREM SAGAR`_
- Skip ``test_to_parquet_with_get`` if no parquet libs available (:pr:`6188`) `Scott Sanderson`_

Documentation
^^^^^^^^^^^^^

- Added documentation for ``distributed.Event`` class (:pr:`6231`) `Nils Braun`_
- Doc write to remote (:pr:`6124`) `Ray Bell`_


.. _v2.16.0 / 2020-05-08:

2.16.0 / 2020-05-08
-------------------

Array
^^^^^

- Fix array general-reduction name (:pr:`6176`) `Nick Evans`_
- Replace ``dim`` with ``shape`` in ``unravel_index`` (:pr:`6155`) `Julia Signell`_
- Moment: handle all elements being masked (:pr:`5339`) `Gabe Joseph`_

Core
^^^^

- Remove Redundant string concatenations in dask code-base (:pr:`6137`) `GALI PREM SAGAR`_
- Upstream compat (:pr:`6159`) `Tom Augspurger`_
- Ensure ``sizeof`` of dict and sequences returns an integer (:pr:`6179`) `James Bourbeau`_
- Estimate python collection sizes with random sampling (:pr:`6154`) `Florian Jetter`_
- Update test upstream (:pr:`6146`) `Tom Augspurger`_
- Skip test for mindeps build (:pr:`6144`) `Tom Augspurger`_
- Switch default multiprocessing context to "spawn" (:pr:`4003`) `Itamar Turner-Trauring`_
- Update manifest to include dask-schema (:pr:`6140`) `Benjamin Zaitlen`_

DataFrame
^^^^^^^^^

- Harden inconsistent-schema handling in pyarrow-based ``read_parquet`` (:pr:`6160`) `Richard (Rick) Zamora`_
- Add compute ``kwargs`` to methods that write data to disk (:pr:`6056`) `Krishan Bhasin`_
- Fix issue where ``unique`` returns an index like result from backends (:pr:`6153`) `GALI PREM SAGAR`_
- Fix internal error in ``map_partitions`` with collections (:pr:`6103`) `Tom Augspurger`_

Documentation
^^^^^^^^^^^^^

- Add phase of computation to index TOC (:pr:`6157`) `Benjamin Zaitlen`_
- Remove unused imports in scheduling script (:pr:`6138`) `James Lamb`_
- Fix indent (:pr:`6147`) `Martin Durant`_
- Add Tom's log config example (:pr:`6143`) `Martin Durant`_


.. _v2.15.0 / 2020-04-24:

2.15.0 / 2020-04-24
-------------------

Array
^^^^^

- Update ``dask.array.from_array`` to warn when passed a Dask collection (:pr:`6122`) `James Bourbeau`_
- Un-numpy like behaviour in ``dask.array.pad`` (:pr:`6042`) `Mark Boer`_
- Add support for ``repeats=0`` in ``da.repeat`` (:pr:`6080`) `James Bourbeau`_

Core
^^^^

- Fix yaml layout for schema (:pr:`6132`) `Benjamin Zaitlen`_
- Configuration Reference (:pr:`6069`) `Benjamin Zaitlen`_
- Add configuration option to turn off task fusion (:pr:`6087`) `Matthew Rocklin`_
- Skip pyarrow on windows (:pr:`6094`) `Tom Augspurger`_
- Set limit to maximum length of fused key (:pr:`6057`) `Lucas Rademaker`_
- Add test against #6062 (:pr:`6072`) `Martin Durant`_
- Bump checkout action to v2 (:pr:`6065`) `James Bourbeau`_

DataFrame
^^^^^^^^^

- Generalize categorical calls to support cudf ``Categorical`` (:pr:`6113`) `GALI PREM SAGAR`_
- Avoid reading ``_metadata`` on every worker (:pr:`6017`) `Richard (Rick) Zamora`_
- Use ``group_split_dispatch`` and ``ignore_index`` in ``apply_concat_apply`` (:pr:`6119`) `Richard (Rick) Zamora`_
- Handle new (dtype) pandas metadata with pyarrow (:pr:`6090`) `Richard (Rick) Zamora`_
- Skip ``test_partition_on_cats_pyarrow`` if pyarrow is not installed (:pr:`6112`) `James Bourbeau`_
- Update DataFrame len to handle columns with the same name (:pr:`6111`) `James Bourbeau`_
- ``ArrowEngine`` bug fixes and test coverage (:pr:`6047`) `Richard (Rick) Zamora`_
- Added mode (:pr:`5958`) `Adam Lewis`_

Documentation
^^^^^^^^^^^^^

- Update "helm install" for helm 3 usage (:pr:`6130`) `JulianWgs`_
- Extend preload documentation (:pr:`6077`) `Matthew Rocklin`_
- Fixed small typo in DataFrame ``map_partitions()`` docstring (:pr:`6115`) `Eugene Huang`_
- Fix typo: "double" should be times, not plus (:pr:`6091`) `David Chudzicki`_
- Fix first line of ``array.random.*`` docs (:pr:`6063`) `Martin Durant`_
- Add section about ``Semaphore`` in distributed (:pr:`6053`) `Florian Jetter`_


.. _v2.14.0 / 2020-04-03:

2.14.0 / 2020-04-03
-------------------

Array
^^^^^

- Added ``np.iscomplexobj`` implementation (:pr:`6045`) `Tom Augspurger`_

Core
^^^^

- Update ``test_rearrange_disk_cleanup_with_exception`` to pass without cloudpickle installed (:pr:`6052`) `James Bourbeau`_
- Fixed flaky ``test-rearrange`` (:pr:`5977`) `Tom Augspurger`_

DataFrame
^^^^^^^^^

- Use ``_meta_nonempty`` for dtype casting in ``stack_partitions`` (:pr:`6061`) `mlondschien`_
- Fix bugs in ``_metadata`` creation and filtering in parquet ``ArrowEngine`` (:pr:`6023`) `Richard (Rick) Zamora`_

Documentation
^^^^^^^^^^^^^

- DOC: Add name caveats (:pr:`6040`) `Tom Augspurger`_


.. _v2.13.0 / 2020-03-25:

2.13.0 / 2020-03-25
-------------------

Array
^^^^^

- Support ``dtype`` and other keyword arguments in ``da.random`` (:pr:`6030`) `Matthew Rocklin`_
- Register support for ``cupy`` sparse ``hstack``/``vstack`` (:pr:`5735`) `Corey J. Nolet`_
- Force ``self.name`` to ``str`` in ``dask.array`` (:pr:`6002`) `Chuanzhu Xu`_

Bag
^^^

- Set ``rename_fused_keys`` to ``None`` by default in ``bag.optimize`` (:pr:`6000`) `Lucas Rademaker`_

Core
^^^^

- Copy dict in ``to_graphviz`` to prevent overwriting (:pr:`5996`) `JulianWgs`_
- Stricter pandas ``xfail`` (:pr:`6024`) `Tom Augspurger`_
- Fix CI failures (:pr:`6013`) `James Bourbeau`_
- Update ``toolz`` to 0.8.2 and use ``tlz`` (:pr:`5997`) `Ryan Grout`_
- Move Windows CI builds to GitHub Actions (:pr:`5862`) `James Bourbeau`_

DataFrame
^^^^^^^^^

- Improve path-related exceptions in ``read_hdf`` (:pr:`6032`) `psimaj`_
- Fix ``dtype`` handling in ``dd.concat`` (:pr:`6006`) `mlondschien`_
- Handle cudf's leftsemi and leftanti joins (:pr:`6025`) `Richard J Zamora`_
- Remove unused ``npartitions`` variable in ``dd.from_pandas`` (:pr:`6019`) `Daniel Saxton`_
- Added shuffle to ``DataFrame.random_split`` (:pr:`5980`) `petiop`_

Documentation
^^^^^^^^^^^^^

- Fix indentation in scheduler-overview docs (:pr:`6022`) `Matthew Rocklin`_
- Update task graphs in optimize docs (:pr:`5928`) `Julia Signell`_
- Optionally get rid of intermediary boxes in visualize, and add more labels (:pr:`5976`) `Julia Signell`_


.. _v2.12.0 / 2020-03-06:

2.12.0 / 2020-03-06
-------------------

Array
^^^^^

- Improve reuse of temporaries with numpy (:pr:`5933`) `Bruce Merry`_
- Make ``map_blocks`` with ``block_info`` produce a ``Blockwise`` (:pr:`5896`) `Bruce Merry`_
- Optimize ``make_blockwise_graph`` (:pr:`5940`) `Bruce Merry`_
- Fix axes ordering in ``da.tensordot`` (:pr:`5975`) `Gil Forsyth`_
- Adds empty mode to ``array.pad`` (:pr:`5931`) `Thomas J. Fan`_

Core
^^^^

- Remove ``toolz.memoize`` dependency in ``dask.utils`` (:pr:`5978`) `Ryan Grout`_
- Close pool leaking subprocess (:pr:`5979`) `Tom Augspurger`_
- Pin ``numpydoc`` to ``0.8.0`` (fix double autoescape) (:pr:`5961`) `Gil Forsyth`_
- Register deterministic tokenization for ``range`` objects (:pr:`5947`) `James Bourbeau`_
- Unpin ``msgpack`` in CI (:pr:`5930`) `JAmes Bourbeau`_
- Ensure dot results are placed in unique files. (:pr:`5937`) `Elliott Sales de Andrade`_
- Add remaining optional dependencies to Travis 3.8 CI build environment (:pr:`5920`) `James Bourbeau`_

DataFrame
^^^^^^^^^

- Skip parquet ``getitem`` optimization for some keys (:pr:`5917`) `Tom Augspurger`_
- Add ``ignore_index`` argument to ``rearrange_by_column`` code path (:pr:`5973`) `Richard J Zamora`_
- Add DataFrame and Series ``memory_usage_per_partition`` methods (:pr:`5971`) `James Bourbeau`_
- ``xfail`` test_describe when using Pandas 0.24.2 (:pr:`5948`) `James Bourbeau`_
- Implement ``dask.dataframe.to_numeric`` (:pr:`5929`) `Julia Signell`_
- Add new error message content when columns are in a different order (:pr:`5927`) `Julia Signell`_
- Use shallow copy for assign operations when possible (:pr:`5740`) `Richard J Zamora`_

Documentation
^^^^^^^^^^^^^

- Changed above to below in ``dask.array.triu`` docs (:pr:`5984`) `Henrik Andersson`_
- Array slicing: fix typo in ``slice_with_int_dask_array`` error message (:pr:`5981`) `Gabe Joseph`_
- Grammar and formatting updates to docstrings (:pr:`5963`) `James Lamb`_
- Update develop doc with conda option (:pr:`5939`) `Ray Bell`_
- Update title of DataFrame extension docs (:pr:`5954`) `James Bourbeau`_
- Fixed typos in documentation (:pr:`5962`) `James Lamb`_
- Add original class or module as a ``kwarg`` on ``_bind_*`` methods (:pr:`5946`) `Julia Signell`_
- Add collect list example (:pr:`5938`) `Ray Bell`_
- Update optimization doc for python 3 (:pr:`5926`) `Julia Signell`_


.. _v2.11.0 / 2020-02-19:

2.11.0 / 2020-02-19
-------------------

Array
^^^^^

- Cache result of ``Array.shape`` (:pr:`5916`) `Bruce Merry`_
- Improve accuracy of ``estimate_graph_size`` for ``rechunk`` (:pr:`5907`) `Bruce Merry`_
- Skip rechunk steps that do not alter chunking (:pr:`5909`) `Bruce Merry`_
- Support ``dtype`` and other ``kwargs`` in ``coarsen`` (:pr:`5903`) `Matthew Rocklin`_
- Push chunk override from ``map_blocks`` into blockwise (:pr:`5895`) `Bruce Merry`_
- Avoid using ``rewrite_blockwise`` for a singleton (:pr:`5890`) `Bruce Merry`_
- Optimize ``slices_from_chunks`` (:pr:`5891`) `Bruce Merry`_
- Avoid unnecessary ``__getitem__`` in ``block()`` when chunks have correct dimensionality (:pr:`5884`) `Thomas Robitaille`_

Bag
^^^

- Add ``include_path`` option for ``dask.bag.read_text`` (:pr:`5836`) `Yifan Gu`_
- Fixes ``ValueError`` in delayed execution of bagged NumPy array (:pr:`5828`) `Surya Avala`_

Core
^^^^

- CI: Pin ``msgpack`` (:pr:`5923`) `Tom Augspurger`_
- Rename ``test_inner`` to ``test_outer`` (:pr:`5922`) `Shiva Raisinghani`_
- ``quote`` should quote dicts too (:pr:`5905`) `Bruce Merry`_
- Register a normalizer for literal (:pr:`5898`) `Bruce Merry`_
- Improve layer name synthesis for non-HLGs (:pr:`5888`) `Bruce Merry`_
- Replace flake8 pre-commit-hook with upstream (:pr:`5892`) `Julia Signell`_
- Call pip as a module to avoid warnings (:pr:`5861`) `Cyril Shcherbin`_
- Close ``ThreadPool`` at exit (:pr:`5852`) `Tom Augspurger`_
- Remove ``dask.dataframe`` import in tokenization code (:pr:`5855`) `James Bourbeau`_

DataFrame
^^^^^^^^^

- Require ``pandas>=0.23`` (:pr:`5883`) `Tom Augspurger`_
- Remove lambda from dataframe aggregation (:pr:`5901`) `Matthew Rocklin`_
- Fix exception chaining in ``dataframe/__init__.py`` (:pr:`5882`) `Ram Rachum`_
- Add support for reductions on empty dataframes (:pr:`5804`) `Shiva Raisinghani`_
- Expose ``sort=`` argument for groupby (:pr:`5801`) `Richard J Zamora`_
- Add ``df.empty`` property (:pr:`5711`) `rockwellw`_
- Use parquet read speed-ups from ``fastparquet.api.paths_to_cats``. (:pr:`5821`) `Igor Gotlibovych`_

Documentation
^^^^^^^^^^^^^

- Deprecate ``doc_wraps`` (:pr:`5912`) `Tom Augspurger`_
- Update array internal design docs for HighLevelGraph era (:pr:`5889`) `Bruce Merry`_
- Move over dashboard connection docs (:pr:`5877`) `Matthew Rocklin`_
- Move prometheus docs from distributed.dask.org (:pr:`5876`) `Matthew Rocklin`_
- Removing duplicated DO block at the end (:pr:`5878`) `K.-Michael Aye`_
- ``map_blocks`` see also (:pr:`5874`) `Tom Augspurger`_
- More derived from (:pr:`5871`) `Julia Signell`_
- Fix typo (:pr:`5866`) `Yetunde Dada`_
- Fix typo in ``cloud.rst`` (:pr:`5860`) `Andrew Thomas`_
- Add note pointing to code of conduct and diversity statement (:pr:`5844`) `Matthew Rocklin`_


.. _v2.10.1 / 2020-01-30:

2.10.1 / 2020-01-30
-------------------

- Fix Pandas 1.0 version comparison (:pr:`5851`) `Tom Augspurger`_
- Fix typo in distributed diagnostics documentation (:pr:`5841`) `Gerrit Holl`_


.. _v2.10.0 / 2020-01-28:

2.10.0 / 2020-01-28
-------------------

- Support for pandas 1.0's new ``BooleanDtype`` and ``StringDtype`` (:pr:`5815`) `Tom Augspurger`_
- Compatibility with pandas 1.0's API breaking changes and deprecations (:pr:`5792`) `Tom Augspurger`_
- Fixed non-deterministic tokenization of some extension-array backed pandas objects (:pr:`5813`) `Tom Augspurger`_
- Fixed handling of dataclass class objects in collections (:pr:`5812`) `Matteo De Wint`_
- Fixed resampling with tz-aware dates when one of the endpoints fell in a non-existent time (:pr:`5807`) `dfonnegra`_
- Delay initial Zarr dataset creation until the computation occurs (:pr:`5797`) `Chris Roat`_
- Use parquet dataset statistics in more cases with the ``pyarrow`` engine (:pr:`5799`) `Richard J Zamora`_
- Fixed exception in ``groupby.std()`` when some of the keys were large integers (:pr:`5737`) `H. Thomson Comer`_


.. _v2.9.2 / 2020-01-16:

2.9.2 / 2020-01-16
------------------

Array
^^^^^

- Unify chunks in ``broadcast_arrays`` (:pr:`5765`) `Matthew Rocklin`_

Core
^^^^

- ``xfail`` CSV encoding tests (:pr:`5791`) `Tom Augspurger`_
- Update order to handle empty dask graph (:pr:`5789`) `James Bourbeau`_
- Redo ``dask.order.order`` (:pr:`5646`) `Erik Welch`_

DataFrame
^^^^^^^^^

- Add transparent compression for on-disk shuffle with ``partd`` (:pr:`5786`) `Christian Wesp`_
- Fix ``repr`` for empty dataframes (:pr:`5781`) `Shiva Raisinghani`_
- Pandas 1.0.0RC0 compat (:pr:`5784`) `Tom Augspurger`_
- Remove buggy assertions (:pr:`5783`) `Tom Augspurger`_
- Pandas 1.0 compat (:pr:`5782`) `Tom Augspurger`_
- Fix bug in pyarrow-based ``read_parquet`` on partitioned datasets (:pr:`5777`) `Richard J Zamora`_
- Compat for pandas 1.0 (:pr:`5779`) `Tom Augspurger`_
- Fix groupby/mean error with with categorical index (:pr:`5776`) `Richard J Zamora`_
- Support empty partitions when performing cumulative aggregation (:pr:`5730`) `Matthew Rocklin`_
- ``set_index`` accepts single-item unnested list (:pr:`5760`) `Wes Roach`_
- Fixed partitioning in set index for ordered ``Categorical`` (:pr:`5715`) `Tom Augspurger`_

Documentation
^^^^^^^^^^^^^

- Note additional use case for ``normalize_token.register`` (:pr:`5766`) `Thomas A Caswell`_
- Update bag ``repartition`` docstring (:pr:`5772`) `Timost`_
- Small typos (:pr:`5771`) `Maarten Breddels`_
- Fix typo in Task Expectations docs (:pr:`5767`) `James Bourbeau`_
- Add docs section on task expectations to graph page (:pr:`5764`) `Devin Petersohn`_


.. _v2.9.1 / 2019-12-27:

2.9.1 / 2019-12-27
------------------

Array
^^^^^

-  Support Array.view with dtype=None (:pr:`5736`) `Anderson Banihirwe`_
-  Add dask.array.nanmedian (:pr:`5684`) `Deepak Cherian`_

Core
^^^^

-  xfail test_temporary_directory on Python 3.8 (:pr:`5734`) `James Bourbeau`_
-  Add support for Python 3.8 (:pr:`5603`) `James Bourbeau`_
-  Use id to dedupe constants in rewrite_blockwise (:pr:`5696`) `Jim Crist`_

DataFrame
^^^^^^^^^

-  Raise error when converting a dask dataframe scalar to a boolean (:pr:`5743`) `James Bourbeau`_
-  Ensure dataframe groupby-variance is greater than zero (:pr:`5728`) `Matthew Rocklin`_
-  Fix DataFrame.__iter__ (:pr:`5719`) `Tom Augspurger`_
-  Support Parquet filters in disjunctive normal form, like PyArrow (:pr:`5656`) `Matteo De Wint`_
-  Auto-detect categorical columns in ArrowEngine-based read_parquet (:pr:`5690`) `Richard J Zamora`_
-  Skip parquet getitem optimization tests if no engine found (:pr:`5697`) `James Bourbeau`_
-  Fix independent optimization of parquet-getitem (:pr:`5613`) `Tom Augspurger`_

Documentation
^^^^^^^^^^^^^

-  Update helm config doc (:pr:`5750`) `Ray Bell`_
-  Link to examples.dask.org in several places (:pr:`5733`) `Tom Augspurger`_
-  Add missing " in performance report example (:pr:`5724`) `James Bourbeau`_
-  Resolve several documentation build warnings (:pr:`5685`) `James Bourbeau`_
-  add info on performance_report (:pr:`5713`) `Benjamin Zaitlen`_
-  Add more docs disclaimers (:pr:`5710`) `Julia Signell`_
-  Fix simple typo: wihout -> without (:pr:`5708`) `Tim Gates`_
-  Update numpydoc dependency (:pr:`5694`) `James Bourbeau`_


.. _v2.9.0 / 2019-12-06:

2.9.0 / 2019-12-06
------------------

Array
^^^^^
- Fix ``da.std`` to work with NumPy arrays (:pr:`5681`) `James Bourbeau`_

Core
^^^^
- Register ``sizeof`` functions for Numba and RMM (:pr:`5668`) `John A Kirkham`_
- Update meeting time (:pr:`5682`) `Tom Augspurger`_

DataFrame
^^^^^^^^^
- Modify ``dd.DataFrame.drop`` to use shallow copy (:pr:`5675`) `Richard J Zamora`_
- Fix bug in ``_get_md_row_groups`` (:pr:`5673`) `Richard J Zamora`_
- Close sqlalchemy engine after querying DB (:pr:`5629`) `Krishan Bhasin`_
- Allow ``dd.map_partitions`` to not enforce meta (:pr:`5660`) `Matthew Rocklin`_
- Generalize ``concat_unindexed_dataframes`` to support cudf-backend (:pr:`5659`) `Richard J Zamora`_
- Add dataframe resample methods (:pr:`5636`) `Benjamin Zaitlen`_
- Compute length of dataframe as length of first column (:pr:`5635`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^
- Doc fixup (:pr:`5665`) `James Bourbeau`_
- Update doc build instructions (:pr:`5640`) `James Bourbeau`_
- Fix ADL link (:pr:`5639`) `Ray Bell`_
- Add documentation build (:pr:`5617`) `James Bourbeau`_


.. _v2.8.1 / 2019-11-22:

2.8.1 / 2019-11-22
------------------

Array
^^^^^
- Use auto rechunking in ``da.rechunk`` if no value given (:pr:`5605`) `Matthew Rocklin`_

Core
^^^^
- Add simple action to activate GH actions (:pr:`5619`) `James Bourbeau`_

DataFrame
^^^^^^^^^
- Fix "file_path_0" bug in ``aggregate_row_groups`` (:pr:`5627`) `Richard J Zamora`_
- Add ``chunksize`` argument to ``read_parquet`` (:pr:`5607`) `Richard J Zamora`_
- Change ``test_repartition_npartitions`` to support arch64 architecture (:pr:`5620`) `ossdev07`_
- Categories lost after groupby + agg (:pr:`5423`) `Oliver Hofkens`_
- Fixed relative path issue with parquet metadata file (:pr:`5608`) `Nuno Gomes Silva`_
- Enable gpu-backed covariance/correlation in dataframes (:pr:`5597`) `Richard J Zamora`_

Documentation
^^^^^^^^^^^^^
- Fix institutional faq and unknown doc warnings (:pr:`5616`) `James Bourbeau`_
- Add doc for some utils (:pr:`5609`) `Tom Augspurger`_
- Removes ``html_extra_path`` (:pr:`5614`) `James Bourbeau`_
- Fixed See Also referencence (:pr:`5612`) `Tom Augspurger`_


.. _v2.8.0 / 2019-11-14:

2.8.0 / 2019-11-14
------------------

Array
^^^^^
-  Implement complete dask.array.tile function (:pr:`5574`) `Bouwe Andela`_
-  Add median along an axis with automatic rechunking (:pr:`5575`) `Matthew Rocklin`_
-  Allow da.asarray to chunk inputs (:pr:`5586`) `Matthew Rocklin`_

Bag
^^^

-  Use key_split in Bag name (:pr:`5571`) `Matthew Rocklin`_

Core
^^^^
-  Switch Doctests to Py3.7 (:pr:`5573`) `Ryan Nazareth`_
-  Relax get_colors test to adapt to new Bokeh release (:pr:`5576`) `Matthew Rocklin`_
-  Add dask.blockwise.fuse_roots optimization (:pr:`5451`) `Matthew Rocklin`_
-  Add sizeof implementation for small dicts (:pr:`5578`) `Matthew Rocklin`_
-  Update fsspec, gcsfs, s3fs (:pr:`5588`) `Tom Augspurger`_

DataFrame
^^^^^^^^^
-  Add dropna argument to groupby (:pr:`5579`) `Richard J Zamora`_
-  Revert "Remove import of dask_cudf, which is now a part of cudf (:pr:`5568`)" (:pr:`5590`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^

-  Add best practice for dask.compute function (:pr:`5583`) `Matthew Rocklin`_
-  Create FUNDING.yml (:pr:`5587`) `Gina Helfrich`_
-  Add screencast for coordination primitives (:pr:`5593`) `Matthew Rocklin`_
-  Move funding to .github repo (:pr:`5589`) `Tom Augspurger`_
-  Update calendar link (:pr:`5569`) `Tom Augspurger`_


.. _v2.7.0 / 2019-11-08:

2.7.0 / 2019-11-08
------------------

This release drops support for Python 3.5

Array
^^^^^

-  Reuse code for assert_eq util method (:pr:`5496`) `Vijayant`_
-  Update da.array to always return a dask array (:pr:`5510`) `James Bourbeau`_
-  Skip transpose on trivial inputs (:pr:`5523`) `Ryan Abernathey`_
-  Avoid NumPy scalar string representation in tokenize (:pr:`5527`) `James Bourbeau`_
-  Remove unnecessary tiledb shape constraint (:pr:`5545`) `Norman Barker`_
-  Removes bytes from sparse array HTML repr (:pr:`5556`) `James Bourbeau`_

Core
^^^^

-  Drop Python 3.5 (:pr:`5528`) `James Bourbeau`_
-  Update the use of fixtures in distributed tests (:pr:`5497`) `Matthew Rocklin`_
-  Changed deprecated bokeh-port to dashboard-address (:pr:`5507`) `darindf`_
-  Avoid updating with identical dicts in ensure_dict (:pr:`5501`) `James Bourbeau`_
-  Test Upstream (:pr:`5516`) `Tom Augspurger`_
-  Accelerate reverse_dict (:pr:`5479`) `Ryan Grout`_
-  Update test_imports.sh (:pr:`5534`) `James Bourbeau`_
-  Support cgroups limits on cpu count in multiprocess and threaded schedulers (:pr:`5499`) `Albert DeFusco`_
-  Update minimum pyarrow version on CI (:pr:`5562`) `James Bourbeau`_
-  Make cloudpickle optional (:pr:`5511`) `Guido Imperiale`_

DataFrame
^^^^^^^^^

-  Add an example of index_col usage (:pr:`3072`) `Bruno Bonfils`_
-  Explicitly use iloc for row indexing (:pr:`5500`) `Krishan Bhasin`_
-  Accept dask arrays on columns assignemnt (:pr:`5224`) `Henrique Ribeiro`-
-  Implement unique and value_counts for SeriesGroupBy (:pr:`5358`) `Scott Sievert`_
-  Add sizeof definition for pyarrow tables and columns (:pr:`5522`) `Richard J Zamora`_
-  Enable row-group task partitioning in pyarrow-based read_parquet (:pr:`5508`) `Richard J Zamora`_
-  Removes npartitions='auto' from dd.merge docstring (:pr:`5531`) `James Bourbeau`_
-  Apply enforce error message shows non-overlapping columns. (:pr:`5530`) `Tom Augspurger`_
-  Optimize meta_nonempty for repetitive dtypes (:pr:`5553`) `Petio Petrov`_
-  Remove import of dask_cudf, which is now a part of cudf (:pr:`5568`) `Mads R. B. Kristensen`_

Documentation
^^^^^^^^^^^^^

-  Make capitalization more consistent in FAQ docs (:pr:`5512`) `Matthew Rocklin`_
-  Add CONTRIBUTING.md (:pr:`5513`) `Jacob Tomlinson`_
-  Document optional dependencies (:pr:`5456`) `Prithvi MK`_
-  Update helm chart docs to reflect new chart repo (:pr:`5539`) `Jacob Tomlinson`_
-  Add Resampler to API docs (:pr:`5551`) `James Bourbeau`_
-  Fix typo in read_sql_table (:pr:`5554`) `Eric Dill`_
-  Add adaptive deployments screencast [skip ci] (:pr:`5566`) `Matthew Rocklin`_


.. _v2.6.0 / 2019-10-15:

2.6.0 / 2019-10-15
------------------

Core
^^^^

- Call ``ensure_dict`` on graphs before entering ``toolz.merge`` (:pr:`5486`) `Matthew Rocklin`_
- Consolidating hash dispatch functions (:pr:`5476`) `Richard J Zamora`_

DataFrame
^^^^^^^^^

- Support Python 3.5 in Parquet code (:pr:`5491`) `Benjamin Zaitlen`_
- Avoid identity check in ``warn_dtype_mismatch`` (:pr:`5489`) `Tom Augspurger`_
- Enable unused groupby tests (:pr:`3480`) `Jörg Dietrich`_
- Remove old parquet and bcolz dataframe optimizations (:pr:`5484`) `Matthew Rocklin`_
- Add getitem optimization for ``read_parquet`` (:pr:`5453`) `Tom Augspurger`_
- Use ``_constructor_sliced`` method to determine Series type (:pr:`5480`) `Richard J Zamora`_
- Fix map(series) for unsorted base series index (:pr:`5459`) `Justin Waugh`_
- Fix ``KeyError`` with Groupby label (:pr:`5467`) `Ryan Nazareth`_

Documentation
^^^^^^^^^^^^^

- Use Zoom meeting instead of appear.in (:pr:`5494`) `Matthew Rocklin`_
- Added curated list of resources (:pr:`5460`) `Javad`_
- Update SSH docs to include ``SSHCluster`` (:pr:`5482`) `Matthew Rocklin`_
- Update "Why Dask?" page (:pr:`5473`) `Matthew Rocklin`_
- Fix typos in docstrings (:pr:`5469`) `garanews`_


.. _v2.5.2 / 2019-10-04:

2.5.2 / 2019-10-04
------------------

Array
^^^^^

-  Correct chunk size logic for asymmetric overlaps (:pr:`5449`) `Ben Jeffery`_
-  Make da.unify_chunks public API (:pr:`5443`) `Matthew Rocklin`_

DataFrame
^^^^^^^^^

-  Fix dask.dataframe.fillna handling of Scalar object (:pr:`5463`) `Zhenqing Li`_

Documentation
^^^^^^^^^^^^^

-  Remove boxes in Spark comparison page (:pr:`5445`) `Matthew Rocklin`_
-  Add latest presentations (:pr:`5446`) `Javad`_
-  Update cloud documentation (:pr:`5444`) `Matthew Rocklin`_


.. _v2.5.0 / 2019-09-27:

2.5.0 / 2019-09-27
------------------

Core
^^^^

-  Add sentinel no_default to get_dependencies task (:pr:`5420`) `James Bourbeau`_
-  Update fsspec version (:pr:`5415`) `Matthew Rocklin`_
-  Remove PY2 checks (:pr:`5400`) `Jim Crist`_

DataFrame
^^^^^^^^^

-  Add option to not check meta in dd.from_delayed (:pr:`5436`) `Christopher J. Wright`_
-  Fix test_timeseries_nulls_in_schema failures with pyarrow master (:pr:`5421`) `Richard J Zamora`_
-  Reduce read_metadata output size in pyarrow/parquet (:pr:`5391`) `Richard J Zamora`_
-  Test numeric edge case for repartition with npartitions. (:pr:`5433`) `amerkel2`_
-  Unxfail pandas-datareader test (:pr:`5430`) `Tom Augspurger`_
-  Add DataFrame.pop implementation (:pr:`5422`) `Matthew Rocklin`_
-  Enable merge/set_index for cudf-based dataframes with cupy ``values`` (:pr:`5322`) `Richard J Zamora`_
-  drop_duplicates support for positional subset parameter (:pr:`5410`) `Wes Roach`_

Documentation
^^^^^^^^^^^^^

-  Add screencasts to array, bag, dataframe, delayed, futures and setup  (:pr:`5429`) (:pr:`5424`) `Matthew Rocklin`_
-  Fix delimeter parsing documentation (:pr:`5428`) `Mahmut Bulut`_
-  Update overview image (:pr:`5404`) `James Bourbeau`_


.. _v2.4.0 / 2019-09-13:

2.4.0 / 2019-09-13
------------------

Array
^^^^^

- Adds explicit ``h5py.File`` mode (:pr:`5390`) `James Bourbeau`_
- Provides method to compute unknown array chunks sizes (:pr:`5312`) `Scott Sievert`_
- Ignore runtime warning in Array ``compute_meta`` (:pr:`5356`) `estebanag`_
- Add ``_meta`` to ``Array.__dask_postpersist__`` (:pr:`5353`) `Benoit Bovy`_
- Fixup ``da.asarray`` and ``da.asanyarray`` for datetime64 dtype and xarray objects (:pr:`5334`) `Stephan Hoyer`_
- Add shape implementation (:pr:`5293`) `Tom Augspurger`_
- Add chunktype to array text repr (:pr:`5289`) `James Bourbeau`_
- Array.random.choice: handle array-like non-arrays (:pr:`5283`) `Gabe Joseph`_

Core
^^^^

- Remove deprecated code (:pr:`5401`) `Jim Crist`_
- Fix ``funcname`` when vectorized func has no ``__name__`` (:pr:`5399`) `James Bourbeau`_
- Truncate ``funcname`` to avoid long key names (:pr:`5383`) `Matthew Rocklin`_
- Add support for ``numpy.vectorize`` in ``funcname`` (:pr:`5396`) `James Bourbeau`_
- Fixed HDFS upstream test (:pr:`5395`) `Tom Augspurger`_
- Support numbers and None in ``parse_bytes``/``timedelta`` (:pr:`5384`) `Matthew Rocklin`_
- Fix tokenizing of subindexes on memmapped numpy arrays (:pr:`5351`) `Henry Pinkard`_
- Upstream fixups (:pr:`5300`) `Tom Augspurger`_

DataFrame
^^^^^^^^^

- Allow pandas to cast type of statistics (:pr:`5402`) `Richard J Zamora`_
- Preserve index dtype after applying ``dd.pivot_table`` (:pr:`5385`) `therhaag`_
- Implement explode for Series and DataFrame (:pr:`5381`) `Arpit Solanki`_
- ``set_index`` on categorical fails with less categories than partitions (:pr:`5354`) `Oliver Hofkens`_
- Support output to a single CSV file (:pr:`5304`) `Hongjiu Zhang`_
- Add ``groupby().transform()`` (:pr:`5327`) `Oliver Hofkens`_
- Adding filter kwarg to pyarrow dataset call (:pr:`5348`) `Richard J Zamora`_
- Implement and check compression defaults for parquet (:pr:`5335`) `Sarah Bird`_
- Pass sqlalchemy params to delayed objects (:pr:`5332`) `Arpit Solanki`_
- Fixing schema handling in arrow-parquet (:pr:`5307`) `Richard J Zamora`_
- Add support for DF and Series ``groupby().idxmin/max()`` (:pr:`5273`) `Oliver Hofkens`_
- Add correlation calculation and add test (:pr:`5296`) `Benjamin Zaitlen`_

Documentation
^^^^^^^^^^^^^

- Numpy docstring standard has moved (:pr:`5405`) `Wes Roach`_
- Reference correct NumPy array name (:pr:`5403`) `Wes Roach`_
- Minor edits to Array chunk documentation (:pr:`5372`) `Scott Sievert`_
- Add methods to API docs (:pr:`5387`) `Tom Augspurger`_
- Add namespacing to configuration example (:pr:`5374`) `Matthew Rocklin`_
- Add get_task_stream and profile to the diagnostics page (:pr:`5375`) `Matthew Rocklin`_
- Add best practice to load data with Dask (:pr:`5369`) `Matthew Rocklin`_
- Update ``institutional-faq.rst`` (:pr:`5345`) `DomHudson`_
- Add threads and processes note to the best practices (:pr:`5340`) `Matthew Rocklin`_
- Update cuDF links (:pr:`5328`) `James Bourbeau`_
- Fixed small typo with parentheses placement (:pr:`5311`) `Eugene Huang`_
- Update link in reshape docstring (:pr:`5297`) `James Bourbeau`_


.. _v2.3.0 / 2019-08-16:

2.3.0 / 2019-08-16
------------------

Array
^^^^^

- Raise exception when ``from_array`` is given a dask array (:pr:`5280`) `David Hoese`_
- Avoid adjusting gufunc's meta dtype twice (:pr:`5274`) `Peter Andreas Entschev`_
- Add ``meta=`` keyword to map_blocks and add test with sparse (:pr:`5269`) `Matthew Rocklin`_
- Add rollaxis and moveaxis (:pr:`4822`) `Tobias de Jong`_
- Always increment old chunk index (:pr:`5256`) `James Bourbeau`_
- Shuffle dask array (:pr:`3901`) `Tom Augspurger`_
- Fix ordering when indexing a dask array with a bool dask array (:pr:`5151`) `James Bourbeau`_

Bag
^^^

- Add workaround for memory leaks in bag generators (:pr:`5208`) `Marco Neumann`_

Core
^^^^

- Set strict xfail option (:pr:`5220`) `James Bourbeau`_
- test-upstream (:pr:`5267`) `Tom Augspurger`_
- Fixed HDFS CI failure (:pr:`5234`) `Tom Augspurger`_
- Error nicely if no file size inferred (:pr:`5231`) `Jim Crist`_
- A few changes to ``config.set`` (:pr:`5226`) `Jim Crist`_
- Fixup black string normalization (:pr:`5227`) `Jim Crist`_
- Pin NumPy in windows tests (:pr:`5228`) `Jim Crist`_
- Ensure parquet tests are skipped if fastparquet and pyarrow not installed (:pr:`5217`) `James Bourbeau`_
- Add fsspec to readthedocs (:pr:`5207`) `Matthew Rocklin`_
- Bump NumPy and Pandas to 1.17 and 0.25 in CI test (:pr:`5179`) `John A Kirkham`_

DataFrame
^^^^^^^^^

- Fix ``DataFrame.query`` docstring (incorrect numexpr API) (:pr:`5271`) `Doug Davis`_
- Parquet metadata-handling improvements (:pr:`5218`) `Richard J Zamora`_
- Improve messaging around sorted parquet columns for index (:pr:`5265`) `Martin Durant`_
- Add ``rearrange_by_divisions`` and ``set_index`` support for cudf (:pr:`5205`) `Richard J Zamora`_
- Fix ``groupby.std()`` with integer colum names (:pr:`5096`) `Nicolas Hug`_
- Add ``Series.__iter__`` (:pr:`5071`) `Blane`_
- Generalize ``hash_pandas_object`` to work for non-pandas backends (:pr:`5184`) `GALI PREM SAGAR`_
- Add rolling cov (:pr:`5154`) `Ivars Geidans`_
- Add columns argument in drop function (:pr:`5223`) `Henrique Ribeiro`_

Documentation
^^^^^^^^^^^^^

- Update institutional FAQ doc (:pr:`5277`) `Matthew Rocklin`_
- Add draft of institutional FAQ (:pr:`5214`) `Matthew Rocklin`_
- Make boxes for dask-spark page (:pr:`5249`) `Martin Durant`_
- Add motivation for shuffle docs (:pr:`5213`) `Matthew Rocklin`_
- Fix links and API entries for best-practices (:pr:`5246`) `Martin Durant`_
- Remove "bytes" (internal data ingestion) doc page (:pr:`5242`) `Martin Durant`_
- Redirect from our local distributed page to distributed.dask.org (:pr:`5248`) `Matthew Rocklin`_
- Cleanup API page (:pr:`5247`) `Matthew Rocklin`_
- Remove excess endlines from install docs (:pr:`5243`) `Matthew Rocklin`_
- Remove item list in phases of computation doc (:pr:`5245`) `Martin Durant`_
- Remove custom graphs from the TOC sidebar (:pr:`5241`) `Matthew Rocklin`_
- Remove experimental status of custom collections (:pr:`5236`) `James Bourbeau`_
- Adds table of contents to Why Dask? (:pr:`5244`) `James Bourbeau`_
- Moves bag overview to top-level bag page (:pr:`5240`) `James Bourbeau`_
- Remove use-cases in favor of stories.dask.org (:pr:`5238`) `Matthew Rocklin`_
- Removes redundant TOC information in index.rst (:pr:`5235`) `James Bourbeau`_
- Elevate dashboard in distributed diagnostics documentation (:pr:`5239`) `Martin Durant`_
- Updates "add" layer in HLG docs example (:pr:`5237`) `James Bourbeau`_
- Update GUFunc documentation (:pr:`5232`) `Matthew Rocklin`_


.. _v2.2.0 / 2019-08-01:

2.2.0 / 2019-08-01
------------------

Array
^^^^^

-  Use da.from_array(..., asarray=False) if input follows NEP-18 (:pr:`5074`) `Matthew Rocklin`_
-  Add missing attributes to from_array documentation (:pr:`5108`) `Peter Andreas Entschev`_
-  Fix meta computation for some reduction functions (:pr:`5035`) `Peter Andreas Entschev`_
-  Raise informative error in to_zarr if unknown chunks (:pr:`5148`) `James Bourbeau`_
-  Remove invalid pad tests (:pr:`5122`) `Tom Augspurger`_
-  Ignore NumPy warnings in compute_meta (:pr:`5103`) `Peter Andreas Entschev`_
-  Fix kurtosis calc for single dimension input array (:pr:`5177`) `@andrethrill`_
-  Support Numpy 1.17 in tests (:pr:`5192`) `Matthew Rocklin`_

Bag
^^^

-  Supply pool to bag test to resolve intermittent failure (:pr:`5172`) `Tom Augspurger`_

Core
^^^^

-  Base dask on fsspec (:pr:`5064`) (:pr:`5121`) `Martin Durant`_
-  Various upstream compatibility fixes (:pr:`5056`) `Tom Augspurger`_
-  Make distributed tests optional again. (:pr:`5128`) `Elliott Sales de Andrade`_
-  Fix HDFS in dask (:pr:`5130`) `Martin Durant`_
-  Ignore some more invalid value warnings. (:pr:`5140`) `Elliott Sales de Andrade`_

DataFrame
^^^^^^^^^

-  Fix pd.MultiIndex size estimate (:pr:`5066`) `Brett Naul`_
-  Generalizing has_known_categories (:pr:`5090`) `GALI PREM SAGAR`_
-  Refactor Parquet engine (:pr:`4995`) `Richard J Zamora`_
-  Add divide method to series and dataframe (:pr:`5094`) `msbrown47`_
-  fix flaky partd test (:pr:`5111`) `Tom Augspurger`_
-  Adjust is_dataframe_like to adjust for value_counts change (:pr:`5143`) `Tom Augspurger`_
-  Generalize rolling windows to support non-Pandas dataframes (:pr:`5149`) `Nick Becker`_
-  Avoid unnecessary aggregation in pivot_table (:pr:`5173`) `Daniel Saxton`_
-  Add column names to apply_and_enforce error message (:pr:`5180`) `Matthew Rocklin`_
-  Add schema keyword argument to to_parquet (:pr:`5150`) `Sarah Bird`_
-  Remove recursion error in accessors (:pr:`5182`) `Jim Crist`_
-  Allow fastparquet to handle gather_statistics=False for file lists (:pr:`5157`) `Richard J Zamora`_

Documentation
^^^^^^^^^^^^^

-  Adds NumFOCUS badge to the README (:pr:`5086`) `James Bourbeau`_
-  Update developer docs [ci skip] (:pr:`5093`) `Jim Crist`_
-  Document DataFrame.set_index computataion behavior `Natalya Rapstine`_
-  Use pip install . instead of calling setup.py (:pr:`5139`) `Matthias Bussonier`_
-  Close user survey (:pr:`5147`) `Tom Augspurger`_
-  Fix Google Calendar meeting link (:pr:`5155`) `Loïc Estève`_
-  Add docker image customization example (:pr:`5171`) `James Bourbeau`_
-  Update remote-data-services after fsspec (:pr:`5170`) `Martin Durant`_
-  Fix typo in spark.rst (:pr:`5164`) `Xavier Holt`_
-  Update setup/python docs for async/await API (:pr:`5163`) `Matthew Rocklin`_
-  Update Local Storage HPC documentation (:pr:`5165`) `Matthew Rocklin`_



.. _v2.1.0 / 2019-07-08:

2.1.0 / 2019-07-08
------------------

Array
^^^^^

- Add ``recompute=`` keyword to ``svd_compressed`` for lower-memory use (:pr:`5041`) `Matthew Rocklin`_
- Change ``__array_function__`` implementation for backwards compatibility (:pr:`5043`) `Ralf Gommers`_
- Added ``dtype`` and ``shape`` kwargs to ``apply_along_axis`` (:pr:`3742`) `Davis Bennett`_
- Fix reduction with empty tuple axis (:pr:`5025`) `Peter Andreas Entschev`_
- Drop size 0 arrays in ``stack`` (:pr:`4978`) `John A Kirkham`_

Core
^^^^

- Removes index keyword from pandas ``to_parquet`` call (:pr:`5075`) `James Bourbeau`_
- Fixes upstream dev CI build installation (:pr:`5072`) `James Bourbeau`_
- Ensure scalar arrays are not rendered to SVG (:pr:`5058`) `Willi Rath`_
- Environment creation overhaul (:pr:`5038`) `Tom Augspurger`_
- s3fs, moto compatibility (:pr:`5033`) `Tom Augspurger`_
- pytest 5.0 compat (:pr:`5027`) `Tom Augspurger`_

DataFrame
^^^^^^^^^

- Fix ``compute_meta`` recursion in blockwise (:pr:`5048`) `Peter Andreas Entschev`_
- Remove hard dependency on pandas in ``get_dummies`` (:pr:`5057`) `GALI PREM SAGAR`_
- Check dtypes unchanged when using ``DataFrame.assign`` (:pr:`5047`) `asmith26`_
- Fix cumulative functions on tables with more than 1 partition (:pr:`5034`) `tshatrov`_
- Handle non-divisible sizes in repartition (:pr:`5013`) `George Sakkis`_
- Handles timestamp and ``preserve_index`` changes in pyarrow (:pr:`5018`) `Richard J Zamora`_
- Fix undefined ``meta`` for ``str.split(expand=False)`` (:pr:`5022`) `Brett Naul`_
- Removed checks used for debugging ``merge_asof`` (:pr:`5011`) `Cody Johnson`_
- Don't use type when getting accessor in dataframes (:pr:`4992`) `Matthew Rocklin`_
- Add ``melt`` as a method of Dask DataFrame (:pr:`4984`) `Dustin Tindall`_
- Adds path-like support to ``to_hdf`` (:pr:`5003`) `James Bourbeau`_

Documentation
^^^^^^^^^^^^^

- Point to latest K8s setup article in JupyterHub docs (:pr:`5065`) `Sean McKenna`_
- Changes vizualize to visualize (:pr:`5061`) `David Brochart`_
- Fix ``from_sequence`` typo in delayed best practices (:pr:`5045`) `James Bourbeau`_
- Add user survey link to docs (:pr:`5026`) `James Bourbeau`_
- Fixes typo in optimization docs (:pr:`5015`) `James Bourbeau`_
- Update community meeting information (:pr:`5006`) `Tom Augspurger`_


.. _v2.0.0 / 2019-06-25:

2.0.0 / 2019-06-25
------------------

Array
^^^^^

-  Support automatic chunking in da.indices (:pr:`4981`) `James Bourbeau`_
-  Err if there are no arrays to stack (:pr:`4975`) `John A Kirkham`_
-  Asymmetrical Array Overlap (:pr:`4863`) `Michael Eaton`_
-  Dispatch concatenate where possible within dask array (:pr:`4669`) `Hameer Abbasi`_
-  Fix tokenization of memmapped numpy arrays on different part of same file (:pr:`4931`) `Henry Pinkard`_
-  Preserve NumPy condition in da.asarray to preserve output shape (:pr:`4945`) `Alistair Miles`_
-  Expand foo_like_safe usage (:pr:`4946`) `Peter Andreas Entschev`_
-  Defer order/casting einsum parameters to NumPy implementation (:pr:`4914`) `Peter Andreas Entschev`_
-  Remove numpy warning in moment calculation (:pr:`4921`) `Matthew Rocklin`_
-  Fix meta_from_array to support Xarray test suite (:pr:`4938`) `Matthew Rocklin`_
-  Cache chunk boundaries for integer slicing (:pr:`4923`) `Bruce Merry`_
-  Drop size 0 arrays in concatenate (:pr:`4167`) `John A Kirkham`_
-  Raise ValueError if concatenate is given no arrays (:pr:`4927`) `John A Kirkham`_
-  Promote types in `concatenate` using `_meta` (:pr:`4925`) `John A Kirkham`_
-  Add chunk type to html repr in Dask array (:pr:`4895`) `Matthew Rocklin`_
-  Add Dask Array._meta attribute (:pr:`4543`) `Peter Andreas Entschev`_
    -  Fix _meta slicing of flexible types (:pr:`4912`) `Peter Andreas Entschev`_
    -  Minor meta construction cleanup in concatenate (:pr:`4937`) `Peter Andreas Entschev`_
    -  Further relax Array meta checks for Xarray (:pr:`4944`) `Matthew Rocklin`_
    -  Support meta= keyword in da.from_delayed (:pr:`4972`) `Matthew Rocklin`_
    -  Concatenate meta along axis (:pr:`4977`) `John A Kirkham`_
    -  Use meta in stack (:pr:`4976`) `John A Kirkham`_
    -  Move blockwise_meta to more general compute_meta function (:pr:`4954`) `Matthew Rocklin`_
-  Alias .partitions to .blocks attribute of dask arrays (:pr:`4853`) `Genevieve Buckley`_
-  Drop outdated `numpy_compat` functions (:pr:`4850`) `John A Kirkham`_
-  Allow da.eye to support arbitrary chunking sizes with chunks='auto'  (:pr:`4834`) `Anderson Banihirwe`_
-  Fix CI warnings in dask.array tests (:pr:`4805`) `Tom Augspurger`_
-  Make map_blocks work with drop_axis + block_info (:pr:`4831`) `Bruce Merry`_
-  Add SVG image and table in Array._repr_html_ (:pr:`4794`) `Matthew Rocklin`_
-  ufunc: avoid __array_wrap__ in favor of __array_function__ (:pr:`4708`) `Peter Andreas Entschev`_
-  Ensure trivial padding returns the original array (:pr:`4990`) `John A Kirkham`_
-  Test ``da.block`` with 0-size arrays (:pr:`4991`) `John A Kirkham`_


Core
^^^^

-  **Drop Python 2.7** (:pr:`4919`) `Jim Crist`_
-  Quiet dependency installs in CI (:pr:`4960`) `Tom Augspurger`_
-  Raise on warnings in tests (:pr:`4916`) `Tom Augspurger`_
-  Add a diagnostics extra to setup.py (includes bokeh) (:pr:`4924`) `John A Kirkham`_
-  Add newline delimter keyword to OpenFile (:pr:`4935`) `btw08`_
-  Overload HighLevelGraphs values method (:pr:`4918`) `James Bourbeau`_
-  Add __await__ method to Dask collections (:pr:`4901`) `Matthew Rocklin`_
-  Also ignore AttributeErrors which may occur if snappy (not python-snappy) is installed (:pr:`4908`) `Mark Bell`_
-  Canonicalize key names in config.rename (:pr:`4903`) `Ian Bolliger`_
-  Bump minimum partd to 0.3.10 (:pr:`4890`) `Tom Augspurger`_
-  Catch async def SyntaxError (:pr:`4836`) `James Bourbeau`_
-  catch IOError in ensure_file (:pr:`4806`) `Justin Poehnelt`_
-  Cleanup CI warnings (:pr:`4798`) `Tom Augspurger`_
-  Move distributed's parse and format functions to dask.utils (:pr:`4793`) `Matthew Rocklin`_
-  Apply black formatting (:pr:`4983`) `James Bourbeau`_
-  Package license file in wheels (:pr:`4988`) `John A Kirkham`_


DataFrame
^^^^^^^^^

-  Add an optional partition_size parameter to repartition (:pr:`4416`) `George Sakkis`_
-  merge_asof and prefix_reduction (:pr:`4877`) `Cody Johnson`_
-  Allow dataframes to be indexed by dask arrays (:pr:`4882`) `Endre Mark Borza`_
-  Avoid deprecated message parameter in pytest.raises (:pr:`4962`) `James Bourbeau`_
-  Update test_to_records to test with lengths argument(:pr:`4515`) `asmith26`_
-  Remove pandas pinning in Dataframe accessors (:pr:`4955`) `Matthew Rocklin`_
-  Fix correlation of series with same names (:pr:`4934`) `Philipp S. Sommer`_
-  Map Dask Series to Dask Series (:pr:`4872`) `Justin Waugh`_
-  Warn in dd.merge on dtype warning (:pr:`4917`) `mcsoini`_
-  Add groupby Covariance/Correlation (:pr:`4889`) `Benjamin Zaitlen`_
-  keep index name with to_datetime (:pr:`4905`) `Ian Bolliger`_
-  Add Parallel variance computation for dataframes (:pr:`4865`) `Ksenia Bobrova`_
-  Add divmod implementation to arrays and dataframes (:pr:`4884`) `Henrique Ribeiro`_
-  Add documentation for dataframe reshape methods (:pr:`4896`) `tpanza`_
-  Avoid use of pandas.compat (:pr:`4881`) `Tom Augspurger`_
-  Added accessor registration for Series, DataFrame, and Index (:pr:`4829`) `Tom Augspurger`_
-  Add read_function keyword to read_json (:pr:`4810`) `Richard J Zamora`_
-  Provide full type name in check_meta (:pr:`4819`) `Matthew Rocklin`_
-  Correctly estimate bytes per row in read_sql_table (:pr:`4807`) `Lijo Jose`_
-  Adding support of non-numeric data to describe() (:pr:`4791`) `Ksenia Bobrova`_
-  Scalars for extension dtypes. (:pr:`4459`) `Tom Augspurger`_
-  Call head before compute in dd.from_delayed (:pr:`4802`) `Matthew Rocklin`_
-  Add support for rolling operations with larger window that partition size in DataFrames with Time-based index (:pr:`4796`) `Jorge Pessoa`_
-  Update groupby-apply doc with warning (:pr:`4800`) `Tom Augspurger`_
-  Change groupby-ness tests in `_maybe_slice` (:pr:`4786`) `Benjamin Zaitlen`_
-  Add master best practices document (:pr:`4745`) `Matthew Rocklin`_
-  Add document for how Dask works with GPUs (:pr:`4792`) `Matthew Rocklin`_
-  Add cli API docs (:pr:`4788`) `James Bourbeau`_
-  Ensure concat output has coherent dtypes (:pr:`4692`) `Guillaume Lemaitre`_
-  Fixes pandas_datareader dependencies installation (:pr:`4989`) `James Bourbeau`_
-  Accept pathlib.Path as pattern in read_hdf (:pr:`3335`) `Jörg Dietrich`_


Documentation
^^^^^^^^^^^^^

-  Move CLI API docs to relavant pages (:pr:`4980`) `James Bourbeau`_
-  Add to_datetime function to dataframe API docs `Matthew Rocklin`_
-  Add documentation entry for dask.array.ma.average (:pr:`4970`) `Bouwe Andela`_
-  Add bag.read_avro to bag API docs (:pr:`4969`) `James Bourbeau`_
-  Fix typo (:pr:`4968`) `mbarkhau`_
-  Docs: Drop support for Python 2.7 (:pr:`4932`) `Hugo`_
-  Remove requirement to modify changelog (:pr:`4915`) `Matthew Rocklin`_
-  Add documentation about meta column order (:pr:`4887`) `Tom Augspurger`_
-  Add documentation note in DataFrame.shift (:pr:`4886`) `Tom Augspurger`_
-  Docs: Fix typo (:pr:`4868`) `Paweł Kordek`_
-  Put do/don't into boxes for delayed best practice docs (:pr:`3821`) `Martin Durant`_
-  Doc fixups (:pr:`2528`) `Tom Augspurger`_
-  Add quansight to paid support doc section (:pr:`4838`) `Martin Durant`_
-  Add document for custom startup (:pr:`4833`) `Matthew Rocklin`_
-  Allow `utils.derive_from` to accept functions, apply across array (:pr:`4804`) `Martin Durant`_
-  Add "Avoid Large Partitions" section to best practices (:pr:`4808`) `Matthew Rocklin`_
-  Update URL for joblib to new website hosting their doc (:pr:`4816`) `Christian Hudon`_

.. _v1.2.2 / 2019-05-08:

1.2.2 / 2019-05-08
------------------

Array
^^^^^

- Clarify regions kwarg to array.store (:pr:`4759`) `Martin Durant`_
- Add dtype= parameter to da.random.randint (:pr:`4753`) `Matthew Rocklin`_
- Use "row major" rather than "C order" in docstring (:pr:`4452`) `@asmith26`_
- Normalize Xarray datasets to Dask arrays (:pr:`4756`) `Matthew Rocklin`_
- Remove normed keyword in da.histogram (:pr:`4755`) `Matthew Rocklin`_

Bag
^^^

- Add key argument to Bag.distinct (:pr:`4423`) `Daniel Severo`_

Core
^^^^

- Add core dask config file (:pr:`4774`) `Matthew Rocklin`_
- Add core dask config file to MANIFEST.in (:pr:`4780`) `James Bourbeau`_
- Enabling glob with HTTP file-system (:pr:`3926`) `Martin Durant`_
- HTTPFile.seek with whence=1 (:pr:`4751`) `Martin Durant`_
- Remove config key normalization (:pr:`4742`) `Jim Crist`_

DataFrame
^^^^^^^^^

- Remove explicit references to Pandas in dask.dataframe.groupby (:pr:`4778`) `Matthew Rocklin`_
- Add support for group_keys kwarg in DataFrame.groupby() (:pr:`4771`) `Brian Chu`_
- Describe doc (:pr:`4762`) `Martin Durant`_
- Remove explicit pandas check in cumulative aggregations (:pr:`4765`) `Nick Becker`_
- Added meta for read_json and test (:pr:`4588`) `Abhinav Ralhan`_
- Add test for dtype casting (:pr:`4760`) `Martin Durant`_
- Document alignment in map_partitions (:pr:`4757`) `Jim Crist`_
- Implement Series.str.split(expand=True) (:pr:`4744`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^

- Tweaks to develop.rst from trying to run tests (:pr:`4772`) `Christian Hudon`_
- Add document describing phases of computation (:pr:`4766`) `Matthew Rocklin`_
- Point users to Dask-Yarn from spark documentation (:pr:`4770`) `Matthew Rocklin`_
- Update images in delayed doc to remove labels (:pr:`4768`) `Martin Durant`_
- Explain intermediate storage for dask arrays (:pr:`4025`) `John A Kirkham`_
- Specify bash code-block in array best practices (:pr:`4764`) `James Bourbeau`_
- Add array best practices doc (:pr:`4705`) `Matthew Rocklin`_
- Update optimization docs now that cull is not automatic (:pr:`4752`) `Matthew Rocklin`_


.. _v1.2.1 / 2019-04-29:

1.2.1 / 2019-04-29
------------------

Array
^^^^^

-  Fix map_blocks with block_info and broadcasting (:pr:`4737`) `Bruce Merry`_
-  Make 'minlength' keyword argument optional in da.bincount (:pr:`4684`) `Genevieve Buckley`_
-  Add support for map_blocks with no array arguments (:pr:`4713`) `Bruce Merry`_
-  Add dask.array.trace (:pr:`4717`) `Danilo Horta`_
-  Add sizeof support for cupy.ndarray (:pr:`4715`) `Peter Andreas Entschev`_
-  Add name kwarg to from_zarr (:pr:`4663`) `Michael Eaton`_
-  Add chunks='auto' to from_array (:pr:`4704`) `Matthew Rocklin`_
-  Raise TypeError if dask array is given as shape for da.ones, zeros, empty or full (:pr:`4707`) `Genevieve Buckley`_
-  Add TileDB backend (:pr:`4679`) `Isaiah Norton`_

Core
^^^^

-  Delay long list arguments (:pr:`4735`) `Matthew Rocklin`_
-  Bump to numpy >= 1.13, pandas >= 0.21.0 (:pr:`4720`) `Jim Crist`_
-  Remove file "test" (:pr:`4710`) `James Bourbeau`_
-  Reenable development build, uses upstream libraries (:pr:`4696`) `Peter Andreas Entschev`_
-  Remove assertion in HighLevelGraph constructor (:pr:`4699`) `Matthew Rocklin`_

DataFrame
^^^^^^^^^

-  Change cum-aggregation last-nonnull-value algorithm (:pr:`4736`) `Nick Becker`_
-  Fixup series-groupby-apply (:pr:`4738`) `Jim Crist`_
-  Refactor array.percentile and dataframe.quantile to use t-digest (:pr:`4677`) `Janne Vuorela`_
-  Allow naive concatenation of sorted dataframes (:pr:`4725`) `Matthew Rocklin`_
-  Fix perf issue in dd.Series.isin (:pr:`4727`) `Jim Crist`_
-  Remove hard pandas dependency for melt by using methodcaller (:pr:`4719`) `Nick Becker`_
-  A few dataframe metadata fixes (:pr:`4695`) `Jim Crist`_
-  Add Dataframe.replace (:pr:`4714`) `Matthew Rocklin`_
-  Add 'threshold' parameter to pd.DataFrame.dropna (:pr:`4625`) `Nathan Matare`_

Documentation
^^^^^^^^^^^^^

-   Add warning about derived docstrings early in the docstring (:pr:`4716`) `Matthew Rocklin`_
-   Create dataframe best practices doc (:pr:`4703`) `Matthew Rocklin`_
-   Uncomment dask_sphinx_theme (:pr:`4728`) `James Bourbeau`_
-   Fix minor typo fix in a Queue/fire_and_forget example (:pr:`4709`) `Matthew Rocklin`_
-   Update from_pandas docstring to match signature (:pr:`4698`) `James Bourbeau`_

.. _v1.2.0 / 2019-04-12:

1.2.0 / 2019-04-12
------------------

Array
^^^^^

-  Fixed mean() and moment() on sparse arrays (:pr:`4525`) `Peter Andreas Entschev`_
-  Add test for NEP-18. (:pr:`4675`) `Hameer Abbasi`_
-  Allow None to say "no chunking" in normalize_chunks (:pr:`4656`) `Matthew Rocklin`_
-  Fix limit value in auto_chunks (:pr:`4645`) `Matthew Rocklin`_

Core
^^^^

-  Updated diagnostic bokeh test for compatibility with bokeh>=1.1.0 (:pr:`4680`) `Philipp Rudiger`_
-  Adjusts codecov's target/threshold, disable patch (:pr:`4671`) `Peter Andreas Entschev`_
-  Always start with empty http buffer, not None (:pr:`4673`) `Martin Durant`_

DataFrame
^^^^^^^^^

-  Propagate index dtype and name when create dask dataframe from array (:pr:`4686`) `Henrique Ribeiro`_
-  Fix ordering of quantiles in describe (:pr:`4647`) `gregrf`_
-  Clean up and document rearrange_column_by_tasks (:pr:`4674`) `Matthew Rocklin`_
-  Mark some parquet tests xfail (:pr:`4667`) `Peter Andreas Entschev`_
-  Fix parquet breakages with arrow 0.13.0 (:pr:`4668`) `Martin Durant`_
-  Allow sample to be False when reading CSV from a remote URL (:pr:`4634`) `Ian Rose`_
-  Fix timezone metadata inference on parquet load (:pr:`4655`) `Martin Durant`_
-  Use is_dataframe/index_like in dd.utils (:pr:`4657`) `Matthew Rocklin`_
-  Add min_count parameter to groupby sum method (:pr:`4648`) `Henrique Ribeiro`_
-  Correct quantile to handle unsorted quantiles (:pr:`4650`) `gregrf`_

Documentation
^^^^^^^^^^^^^

-  Add delayed extra dependencies to install docs (:pr:`4660`) `James Bourbeau`_


.. _v1.1.5 / 2019-03-29:

1.1.5 / 2019-03-29
------------------

Array
^^^^^

-  Ensure that we use the dtype keyword in normalize_chunks (:pr:`4646`) `Matthew Rocklin`_

Core
^^^^

-  Use recursive glob in LocalFileSystem (:pr:`4186`) `Brett Naul`_
-  Avoid YAML deprecation (:pr:`4603`)
-  Fix CI and add set -e (:pr:`4605`) `James Bourbeau`_
-  Support builtin sequence types in dask.visualize (:pr:`4602`)
-  unpack/repack orderedDict (:pr:`4623`) `Justin Poehnelt`_
-  Add da.random.randint to API docs (:pr:`4628`) `James Bourbeau`_
-  Add zarr to CI environment (:pr:`4604`) `James Bourbeau`_
-  Enable codecov (:pr:`4631`) `Peter Andreas Entschev`_

DataFrame
^^^^^^^^^

-  Support setting the index (:pr:`4565`)
-  DataFrame.itertuples accepts index, name kwargs (:pr:`4593`) `Dan O'Donovan`_
-  Support non-Pandas series in dd.Series.unique (:pr:`4599`) `Benjamin Zaitlen`_
-  Replace use of explicit type check with ._is_partition_type predicate (:pr:`4533`)
-  Remove additional pandas warnings in tests (:pr:`4576`)
-  Check object for name/dtype attributes rather than type (:pr:`4606`)
-  Fix comparison against pd.Series (:pr:`4613`) `amerkel2`_
-  Fixing warning from setting categorical codes to floats (:pr:`4624`) `Julia Signell`_
-  Fix renaming on index to_frame method (:pr:`4498`) `Henrique Ribeiro`_
-  Fix divisions when joining two single-partition dataframes (:pr:`4636`) `Justin Waugh`_
-  Warn if partitions overlap in compute_divisions (:pr:`4600`) `Brian Chu`_
-  Give informative meta= warning (:pr:`4637`) `Matthew Rocklin`_
-  Add informative error message to Series.__getitem__ (:pr:`4638`) `Matthew Rocklin`_
-  Add clear exception message when using index or index_col in read_csv (:pr:`4651`) `Álvaro Abella Bascarán`_

Documentation
^^^^^^^^^^^^^

-  Add documentation for custom groupby aggregations (:pr:`4571`)
-  Docs dataframe joins (:pr:`4569`)
-  Specify fork-based contributions  (:pr:`4619`) `James Bourbeau`_
-  correct to_parquet example in docs (:pr:`4641`) `Aaron Fowles`_
-  Update and secure several references (:pr:`4649`) `Søren Fuglede Jørgensen`_


.. _v1.1.4 / 2019-03-08:

1.1.4 / 2019-03-08
------------------

Array
^^^^^

-  Use mask selection in compress (:pr:`4548`) `John A Kirkham`_
-  Use `asarray` in `extract` (:pr:`4549`) `John A Kirkham`_
-  Use correct dtype when test concatenation. (:pr:`4539`) `Elliott Sales de Andrade`_
-  Fix CuPy tests or properly marks as xfail (:pr:`4564`) `Peter Andreas Entschev`_

Core
^^^^

-  Fix local scheduler callback to deal with custom caching (:pr:`4542`) `Yu Feng`_
-  Use parse_bytes in read_bytes(sample=...) (:pr:`4554`) `Matthew Rocklin`_

DataFrame
^^^^^^^^^

-  Fix up groupby-standard deviation again on object dtype keys (:pr:`4541`) `Matthew Rocklin`_
-  TST/CI: Updates for pandas 0.24.1 (:pr:`4551`) `Tom Augspurger`_
-  Add ability to control number of unique elements in timeseries (:pr:`4557`) `Matthew Rocklin`_
-  Add support in read_csv for parameter skiprows for other iterables (:pr:`4560`) `@JulianWgs`_

Documentation
^^^^^^^^^^^^^

-  DataFrame to Array conversion and unknown chunks (:pr:`4516`) `Scott Sievert`_
-  Add docs for random array creation (:pr:`4566`) `Matthew Rocklin`_
-  Fix typo in docstring (:pr:`4572`) `Shyam Saladi`_


.. _v1.1.3 / 2019-03-01:

1.1.3 / 2019-03-01
------------------

Array
^^^^^

-  Modify mean chunk functions to return dicts rather than arrays (:pr:`4513`) `Matthew Rocklin`_
-  Change sparse installation in CI for NumPy/Python2 compatibility (:pr:`4537`) `Matthew Rocklin`_

DataFrame
^^^^^^^^^

-  Make merge dispatchable on pandas/other dataframe types (:pr:`4522`) `Matthew Rocklin`_
-  read_sql_table - datetime index fix and  index type checking (:pr:`4474`) `Joe Corbett`_
-  Use generalized form of index checking (is_index_like) (:pr:`4531`) `Benjamin Zaitlen`_
-  Add tests for groupby reductions with object dtypes (:pr:`4535`) `Matthew Rocklin`_
-  Fixes #4467 : Updates time_series for pandas deprecation (:pr:`4530`) `@HSR05`_

Documentation
^^^^^^^^^^^^^

-  Add missing method to documentation index (:pr:`4528`) `Bart Broere`_


.. _v1.1.2 / 2019-02-25:

1.1.2 / 2019-02-25
------------------

Array
^^^^^

-  Fix another unicode/mixed-type edge case in normalize_array (:pr:`4489`) `Marco Neumann`_
-  Add dask.array.diagonal (:pr:`4431`) `Danilo Horta`_
-  Call asanyarray in unify_chunks (:pr:`4506`) `Jim Crist`_
-  Modify moment chunk functions to return dicts (:pr:`4519`) `Peter Andreas Entschev`_


Bag
^^^

-  Don't inline output keys in dask.bag (:pr:`4464`) `Jim Crist`_
-  Ensure that bag.from_sequence always includes at least one partition (:pr:`4475`) `Anderson Banihirwe`_
-  Implement out_type for bag.fold (:pr:`4502`) `Matthew Rocklin`_
-  Remove map from bag keynames (:pr:`4500`) `Matthew Rocklin`_
-  Avoid itertools.repeat in map_partitions (:pr:`4507`) `Matthew Rocklin`_


DataFrame
^^^^^^^^^

-  Fix relative path parsing on windows when using fastparquet (:pr:`4445`) `Janne Vuorela`_
-  Fix bug in pyarrow and hdfs (:pr:`4453`) (:pr:`4455`) `Michał Jastrzębski`_
-  df getitem with integer slices is not implemented (:pr:`4466`) `Jim Crist`_
-  Replace cudf-specific code with dask-cudf import (:pr:`4470`) `Matthew Rocklin`_
-  Avoid groupby.agg(callable) in groupby-var (:pr:`4482`) `Matthew Rocklin`_
-  Consider uint types as numerical in check_meta (:pr:`4485`) `Marco Neumann`_
-  Fix some typos in groupby comments (:pr:`4494`) `Daniel Saxton`_
-  Add error message around set_index(inplace=True) (:pr:`4501`) `Matthew Rocklin`_
-  meta_nonempty works with categorical index (:pr:`4505`) `Jim Crist`_
-  Add module name to expected meta error message (:pr:`4499`) `Matthew Rocklin`_
-  groupby-nunique works on empty chunk (:pr:`4504`) `Jim Crist`_
-  Propagate index metadata if not specified (:pr:`4509`) `Jim Crist`_

Documentation
^^^^^^^^^^^^^

-  Update docs to use ``from_zarr`` (:pr:`4472`) `John A Kirkham`_
-  DOC: add section of `Using Other S3-Compatible Services` for remote-data-services (:pr:`4405`) `Aploium`_
-  Fix header level of section in changelog (:pr:`4483`) `Bruce Merry`_
-  Add quotes to pip install [skip-ci] (:pr:`4508`) `James Bourbeau`_

Core
^^^^

-  Extend started_cbs AFTER state is initialized (:pr:`4460`) `Marco Neumann`_
-  Fix bug in HTTPFile._fetch_range with headers (:pr:`4479`) (:pr:`4480`) `Ross Petchler`_
-  Repeat optimize_blockwise for diamond fusion (:pr:`4492`) `Matthew Rocklin`_


.. _v1.1.1 / 2019-01-31:

1.1.1 / 2019-01-31
------------------

Array
^^^^^

-  Add support for cupy.einsum (:pr:`4402`) `Johnnie Gray`_
-  Provide byte size in chunks keyword (:pr:`4434`) `Adam Beberg`_
-  Raise more informative error for histogram bins and range (:pr:`4430`) `James Bourbeau`_

DataFrame
^^^^^^^^^

-  Lazily register more cudf functions and move to backends file (:pr:`4396`) `Matthew Rocklin`_
-  Fix ORC tests for pyarrow 0.12.0 (:pr:`4413`) `Jim Crist`_
-  rearrange_by_column: ensure that shuffle arg defaults to 'disk' if it's None in dask.config (:pr:`4414`) `George Sakkis`_
-  Implement filters for _read_pyarrow (:pr:`4415`) `George Sakkis`_
-  Avoid checking against types in is_dataframe_like (:pr:`4418`) `Matthew Rocklin`_
-  Pass username as 'user' when using pyarrow (:pr:`4438`) `Roma Sokolov`_

Delayed
^^^^^^^

-  Fix DelayedAttr return value (:pr:`4440`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^

-  Use SVG for pipeline graphic (:pr:`4406`) `John A Kirkham`_
-  Add doctest-modules to py.test documentation (:pr:`4427`) `Daniel Severo`_

Core
^^^^

-  Work around psutil 5.5.0 not allowing pickling Process objects `Janne Vuorela`_


.. _v1.1.0 / 2019-01-18:

1.1.0 / 2019-01-18
------------------

Array
^^^^^

-  Fix the average function when there is a masked array (:pr:`4236`) `Damien Garaud`_
-  Add allow_unknown_chunksizes to hstack and vstack (:pr:`4287`) `Paul Vecchio`_
-  Fix tensordot for 27+ dimensions (:pr:`4304`) `Johnnie Gray`_
-  Fixed block_info with axes. (:pr:`4301`) `Tom Augspurger`_
-  Use safe_wraps for matmul (:pr:`4346`) `Mark Harfouche`_
-  Use chunks="auto" in array creation routines (:pr:`4354`) `Matthew Rocklin`_
-  Fix np.matmul in dask.array.Array.__array_ufunc__ (:pr:`4363`) `Stephan Hoyer`_
-  COMPAT: Re-enable multifield copy->view change (:pr:`4357`) `Diane Trout`_
-  Calling np.dtype on a delayed object works (:pr:`4387`) `Jim Crist`_
-  Rework normalize_array for numpy data (:pr:`4312`) `Marco Neumann`_

DataFrame
^^^^^^^^^

-  Add fill_value support for series comparisons (:pr:`4250`) `James Bourbeau`_
-  Add schema name in read_sql_table for empty tables (:pr:`4268`) `Mina Farid`_
-  Adjust check for bad chunks in map_blocks (:pr:`4308`) `Tom Augspurger`_
-  Add dask.dataframe.read_fwf (:pr:`4316`) `@slnguyen`_
-  Use atop fusion in dask dataframe (:pr:`4229`) `Matthew Rocklin`_
-  Use parallel_types() in from_pandas (:pr:`4331`) `Matthew Rocklin`_
-  Change DataFrame._repr_data to method (:pr:`4330`) `Matthew Rocklin`_
-  Install pyarrow fastparquet for Appveyor (:pr:`4338`) `Gábor Lipták`_
-  Remove explicit pandas checks and provide cudf lazy registration (:pr:`4359`) `Matthew Rocklin`_
-  Replace isinstance(..., pandas) with is_dataframe_like (:pr:`4375`) `Matthew Rocklin`_
-  ENH: Support 3rd-party ExtensionArrays (:pr:`4379`) `Tom Augspurger`_
-  Pandas 0.24.0 compat (:pr:`4374`) `Tom Augspurger`_

Documentation
^^^^^^^^^^^^^

-  Fix link to 'map_blocks' function in array api docs (:pr:`4258`) `David Hoese`_
-  Add a paragraph on Dask-Yarn in the cloud docs (:pr:`4260`) `Jim Crist`_
-  Copy edit documentation (:pr:`4267`), (:pr:`4263`), (:pr:`4262`), (:pr:`4277`), (:pr:`4271`), (:pr:`4279`), (:pr:`4265`), (:pr:`4295`), (:pr:`4293`), (:pr:`4296`), (:pr:`4302`), (:pr:`4306`), (:pr:`4318`), (:pr:`4314`), (:pr:`4309`), (:pr:`4317`), (:pr:`4326`), (:pr:`4325`), (:pr:`4322`), (:pr:`4332`), (:pr:`4333`), `Miguel Farrajota`_
-  Fix typo in code example (:pr:`4272`) `Daniel Li`_
-  Doc: Update array-api.rst (:pr:`4259`) (:pr:`4282`) `Prabakaran Kumaresshan`_
-  Update hpc doc (:pr:`4266`) `Guillaume Eynard-Bontemps`_
-  Doc: Replace from_avro with read_avro in documents (:pr:`4313`) `Prabakaran Kumaresshan`_
-  Remove reference to "get" scheduler functions in docs (:pr:`4350`) `Matthew Rocklin`_
-  Fix typo in docstring (:pr:`4376`) `Daniel Saxton`_
-  Added documentation for dask.dataframe.merge (:pr:`4382`) `Jendrik Jördening`_

Core
^^^^

-  Avoid recursion in dask.core.get (:pr:`4219`) `Matthew Rocklin`_
-  Remove verbose flag from pytest setup.cfg (:pr:`4281`) `Matthew Rocklin`_
-  Support Pytest 4.0 by specifying marks explicitly (:pr:`4280`) `Takahiro Kojima`_
-  Add High Level Graphs (:pr:`4092`) `Matthew Rocklin`_
-  Fix SerializableLock locked and acquire methods (:pr:`4294`) `Stephan Hoyer`_
-  Pin boto3 to earlier version in tests to avoid moto conflict (:pr:`4276`) `Martin Durant`_
-  Treat None as missing in config when updating (:pr:`4324`) `Matthew Rocklin`_
-  Update Appveyor to Python 3.6 (:pr:`4337`) `Gábor Lipták`_
-  Use parse_bytes more liberally in dask.dataframe/bytes/bag (:pr:`4339`) `Matthew Rocklin`_
-  Add a better error message when cloudpickle is missing (:pr:`4342`) `Mark Harfouche`_
-  Support pool= keyword argument in threaded/multiprocessing get functions (:pr:`4351`) `Matthew Rocklin`_
-  Allow updates from arbitrary Mappings in config.update, not only dicts. (:pr:`4356`) `Stuart Berg`_
-  Move dask/array/top.py code to dask/blockwise.py (:pr:`4348`) `Matthew Rocklin`_
-  Add has_parallel_type (:pr:`4395`) `Matthew Rocklin`_
-  CI: Update Appveyor (:pr:`4381`) `Tom Augspurger`_
-  Ignore non-readable config files (:pr:`4388`) `Jim Crist`_


.. _v1.0.0 / 2018-11-28:

1.0.0 / 2018-11-28
------------------

Array
^^^^^

-  Add nancumsum/nancumprod unit tests (:pr:`4215`) `Guido Imperiale`_

DataFrame
^^^^^^^^^

-  Add index to to_dask_dataframe docstring (:pr:`4232`) `James Bourbeau`_
-  Text and fix when appending categoricals with fastparquet (:pr:`4245`) `Martin Durant`_
-  Don't reread metadata when passing ParquetFile to read_parquet (:pr:`4247`) `Martin Durant`_

Documentation
^^^^^^^^^^^^^

-  Copy edit documentation (:pr:`4222`) (:pr:`4224`) (:pr:`4228`) (:pr:`4231`) (:pr:`4230`) (:pr:`4234`) (:pr:`4235`) (:pr:`4254`) `Miguel Farrajota`_
-  Updated doc for the new scheduler keyword (:pr:`4251`) `@milesial`_


Core
^^^^

-  Avoid a few warnings (:pr:`4223`) `Matthew Rocklin`_
-  Remove dask.store module (:pr:`4221`) `Matthew Rocklin`_
-  Remove AUTHORS.md `Jim Crist`_


.. _v0.20.2 / 2018-11-15:

0.20.2 / 2018-11-15
-------------------

Array
^^^^^

-  Avoid fusing dependencies of atop reductions (:pr:`4207`) `Matthew Rocklin`_

Dataframe
^^^^^^^^^

-  Improve memory footprint for dataframe correlation (:pr:`4193`) `Damien Garaud`_
-  Add empty DataFrame check to boundary_slice (:pr:`4212`) `James Bourbeau`_


Documentation
^^^^^^^^^^^^^

-  Copy edit documentation (:pr:`4197`) (:pr:`4204`) (:pr:`4198`) (:pr:`4199`) (:pr:`4200`) (:pr:`4202`) (:pr:`4209`) `Miguel Farrajota`_
-  Add stats module namespace (:pr:`4206`) `James Bourbeau`_
-  Fix link in dataframe documentation (:pr:`4208`) `James Bourbeau`_


.. _v0.20.1 / 2018-11-09:

0.20.1 / 2018-11-09
-------------------

Array
^^^^^

-  Only allocate the result space in wrapped_pad_func (:pr:`4153`) `John A Kirkham`_
-  Generalize expand_pad_width to expand_pad_value (:pr:`4150`) `John A Kirkham`_
-  Test da.pad with 2D linear_ramp case (:pr:`4162`) `John A Kirkham`_
-  Fix import for broadcast_to. (:pr:`4168`) `samc0de`_
-  Rewrite Dask Array's `pad` to add only new chunks (:pr:`4152`) `John A Kirkham`_
-  Validate index inputs to atop (:pr:`4182`) `Matthew Rocklin`_

Core
^^^^

-  Dask.config set and get normalize underscores and hyphens (:pr:`4143`) `James Bourbeau`_
-  Only subs on core collections, not subclasses (:pr:`4159`) `Matthew Rocklin`_
-  Add block_size=0 option to HTTPFileSystem. (:pr:`4171`) `Martin Durant`_
-  Add traverse support for dataclasses (:pr:`4165`) `Armin Berres`_
-  Avoid optimization on sharedicts without dependencies (:pr:`4181`) `Matthew Rocklin`_
-  Update the pytest version for TravisCI (:pr:`4189`) `Damien Garaud`_
-  Use key_split rather than funcname in visualize names (:pr:`4160`) `Matthew Rocklin`_

Dataframe
^^^^^^^^^

-  Add fix for  DataFrame.__setitem__ for index (:pr:`4151`) `Anderson Banihirwe`_
-  Fix column choice when passing list of files to fastparquet (:pr:`4174`) `Martin Durant`_
-  Pass engine_kwargs from read_sql_table to sqlalchemy (:pr:`4187`) `Damien Garaud`_

Documentation
^^^^^^^^^^^^^

-  Fix documentation in Delayed best practices example that returned an empty list (:pr:`4147`) `Jonathan Fraine`_
-  Copy edit documentation (:pr:`4164`) (:pr:`4175`) (:pr:`4185`) (:pr:`4192`) (:pr:`4191`) (:pr:`4190`) (:pr:`4180`) `Miguel Farrajota`_
-  Fix typo in docstring (:pr:`4183`) `Carlos Valiente`_


.. _v0.20.0 / 2018-10-26:

0.20.0 / 2018-10-26
-------------------

Array
^^^^^

-  Fuse Atop operations (:pr:`3998`), (:pr:`4081`) `Matthew Rocklin`_
-  Support da.asanyarray on dask dataframes (:pr:`4080`) `Matthew Rocklin`_
-  Remove unnecessary endianness check in datetime test (:pr:`4113`) `Elliott Sales de Andrade`_
-  Set name=False in array foo_like functions (:pr:`4116`) `Matthew Rocklin`_
-  Remove dask.array.ghost module (:pr:`4121`) `Matthew Rocklin`_
-  Fix use of getargspec in dask array (:pr:`4125`) `Stephan Hoyer`_
-  Adds dask.array.invert (:pr:`4127`), (:pr:`4131`) `Anderson Banihirwe`_
-  Raise informative error on arg-reduction on unknown chunksize (:pr:`4128`), (:pr:`4135`) `Matthew Rocklin`_
-  Normalize reversed slices in dask array (:pr:`4126`) `Matthew Rocklin`_

Bag
^^^

-  Add bag.to_avro (:pr:`4076`) `Martin Durant`_

Core
^^^^

-  Pull num_workers from config.get (:pr:`4086`), (:pr:`4093`) `James Bourbeau`_
-  Fix invalid escape sequences with raw strings (:pr:`4112`) `Elliott Sales de Andrade`_
-  Raise an error on the use of the get= keyword and set_options (:pr:`4077`) `Matthew Rocklin`_
-  Add import for Azure DataLake storage, and add docs (:pr:`4132`) `Martin Durant`_
-  Avoid collections.Mapping/Sequence (:pr:`4138`)  `Matthew Rocklin`_

Dataframe
^^^^^^^^^

-  Include index keyword in to_dask_dataframe (:pr:`4071`) `Matthew Rocklin`_
-  add support for duplicate column names (:pr:`4087`) `Jan Koch`_
-  Implement min_count for the DataFrame methods sum and prod (:pr:`4090`) `Bart Broere`_
-  Remove pandas warnings in concat (:pr:`4095`) `Matthew Rocklin`_
-  DataFrame.to_csv header option to only output headers in the first chunk (:pr:`3909`) `Rahul Vaidya`_
-  Remove Series.to_parquet (:pr:`4104`) `Justin Dennison`_
-  Avoid warnings and deprecated pandas methods (:pr:`4115`) `Matthew Rocklin`_
-  Swap 'old' and 'previous' when reporting append error (:pr:`4130`) `Martin Durant`_

Documentation
^^^^^^^^^^^^^

-  Copy edit documentation (:pr:`4073`), (:pr:`4074`), (:pr:`4094`), (:pr:`4097`), (:pr:`4107`), (:pr:`4124`), (:pr:`4133`), (:pr:`4139`) `Miguel Farrajota`_
-  Fix typo in code example (:pr:`4089`) `Antonino Ingargiola`_
-  Add pycon 2018 presentation (:pr:`4102`) `Javad`_
-  Quick description for gcsfs (:pr:`4109`) `Martin Durant`_
-  Fixed typo in docstrings of read_sql_table method (:pr:`4114`) `TakaakiFuruse`_
-  Make target directories in redirects if they don't exist (:pr:`4136`) `Matthew Rocklin`_



.. _v0.19.4 / 2018-10-09:

0.19.4 / 2018-10-09
-------------------

Array
^^^^^

-  Implement ``apply_gufunc(..., axes=..., keepdims=...)`` (:pr:`3985`) `Markus Gonser`_

Bag
^^^

-  Fix typo in datasets.make_people (:pr:`4069`) `Matthew Rocklin`_

Dataframe
^^^^^^^^^

-  Added `percentiles` options for `dask.dataframe.describe` method (:pr:`4067`) `Zhenqing Li`_
-  Add DataFrame.partitions accessor similar to Array.blocks (:pr:`4066`) `Matthew Rocklin`_

Core
^^^^

-  Pass get functions and Clients through scheduler keyword (:pr:`4062`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^

-  Fix Typo on hpc example. (missing `=` in kwarg). (:pr:`4068`) `Matthias Bussonier`_
-  Extensive copy-editing: (:pr:`4065`), (:pr:`4064`), (:pr:`4063`) `Miguel Farrajota`_


.. _v0.19.3 / 2018-10-05:

0.19.3 / 2018-10-05
-------------------

Array
^^^^^

-   Make da.RandomState extensible to other modules (:pr:`4041`) `Matthew Rocklin`_
-   Support unknown dims in ravel no-op case (:pr:`4055`) `Jim Crist`_
-   Add basic infrastructure for cupy (:pr:`4019`) `Matthew Rocklin`_
-   Avoid asarray and lock arguments for from_array(getitem) (:pr:`4044`) `Matthew Rocklin`_
-   Move local imports in `corrcoef` to global imports (:pr:`4030`) `John A Kirkham`_
-   Move local `indices` import to global import (:pr:`4029`) `John A Kirkham`_
-   Fix-up Dask Array's fromfunction w.r.t. dtype and kwargs (:pr:`4028`) `John A Kirkham`_
-   Don't use dummy expansion for trim_internal in overlapped (:pr:`3964`) `Mark Harfouche`_
-   Add unravel_index (:pr:`3958`) `John A Kirkham`_

Bag
^^^

-   Sort result in Bag.frequencies (:pr:`4033`) `Matthew Rocklin`_
-   Add support for npartitions=1 edge case in groupby (:pr:`4050`) `James Bourbeau`_
-   Add new random dataset for people (:pr:`4018`) `Matthew Rocklin`_
-   Improve performance of bag.read_text on small files (:pr:`4013`) `Eric Wolak`_
-   Add bag.read_avro (:pr:`4000`) (:pr:`4007`) `Martin Durant`_

Dataframe
^^^^^^^^^

-   Added an ``index`` parameter to :meth:`dask.dataframe.from_dask_array` for creating a dask DataFrame from a dask Array with a given index. (:pr:`3991`) `Tom Augspurger`_
-   Improve sub-classability of dask dataframe (:pr:`4015`) `Matthew Rocklin`_
-   Fix failing hdfs test [test-hdfs] (:pr:`4046`) `Jim Crist`_
-   fuse_subgraphs works without normal fuse (:pr:`4042`) `Jim Crist`_
-   Make path for reading many parquet files without prescan (:pr:`3978`) `Martin Durant`_
-   Index in dd.from_dask_array (:pr:`3991`) `Tom Augspurger`_
-   Making skiprows accept lists (:pr:`3975`) `Julia Signell`_
-   Fail early in fastparquet read for nonexistent column (:pr:`3989`) `Martin Durant`_

Core
^^^^

-   Add support for npartitions=1 edge case in groupby (:pr:`4050`) `James Bourbeau`_
-   Automatically wrap large arguments with dask.delayed in map_blocks/partitions (:pr:`4002`) `Matthew Rocklin`_
-   Fuse linear chains of subgraphs (:pr:`3979`) `Jim Crist`_
-   Make multiprocessing context configurable (:pr:`3763`) `Itamar Turner-Trauring`_

Documentation
^^^^^^^^^^^^^

-   Extensive copy-editing  (:pr:`4049`), (:pr:`4034`),  (:pr:`4031`), (:pr:`4020`), (:pr:`4021`), (:pr:`4022`), (:pr:`4023`), (:pr:`4016`), (:pr:`4017`), (:pr:`4010`), (:pr:`3997`), (:pr:`3996`), `Miguel Farrajota`_
-   Update shuffle method selection docs (:pr:`4048`) `James Bourbeau`_
-   Remove docs/source/examples, point to examples.dask.org (:pr:`4014`) `Matthew Rocklin`_
-   Replace readthedocs links with dask.org (:pr:`4008`) `Matthew Rocklin`_
-   Updates DataFrame.to_hdf docstring for returned values (:pr:`3992`) `James Bourbeau`_


.. _v0.19.2 / 2018-09-17:

0.19.2 / 2018-09-17
-------------------

Array
^^^^^

-  ``apply_gufunc`` implements automatic infer of functions output dtypes (:pr:`3936`) `Markus Gonser`_
-  Fix array histogram range error when array has nans (:pr:`3980`) `James Bourbeau`_
-  Issue 3937 follow up, int type checks. (:pr:`3956`) `Yu Feng`_
-  from_array: add @martindurant's explaining of how hashing is done for an array. (:pr:`3965`) `Mark Harfouche`_
-  Support gradient with coordinate (:pr:`3949`) `Keisuke Fujii`_

Core
^^^^

-  Fix use of has_keyword with partial in Python 2.7 (:pr:`3966`) `Mark Harfouche`_
-  Set pyarrow as default for HDFS (:pr:`3957`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^

-  Use dask_sphinx_theme (:pr:`3963`) `Matthew Rocklin`_
-  Use JupyterLab in Binder links from main page `Matthew Rocklin`_
-  DOC: fixed sphinx syntax (:pr:`3960`) `Tom Augspurger`_


.. _v0.19.1 / 2018-09-06:

0.19.1 / 2018-09-06
-------------------

Array
^^^^^

-  Don't enforce dtype if result has no dtype (:pr:`3928`) `Matthew Rocklin`_
-  Fix NumPy issubtype deprecation warning (:pr:`3939`) `Bruce Merry`_
-  Fix arg reduction tokens to be unique with different arguments (:pr:`3955`) `Tobias de Jong`_
-  Coerce numpy integers to ints in slicing code (:pr:`3944`) `Yu Feng`_
-  Linalg.norm ndim along axis partial fix (:pr:`3933`) `Tobias de Jong`_

Dataframe
^^^^^^^^^

-  Deterministic DataFrame.set_index (:pr:`3867`) `George Sakkis`_
-  Fix divisions in read_parquet when dealing with filters #3831 #3930 (:pr:`3923`) (:pr:`3931`)  `@andrethrill`_
-  Fixing returning type in categorical.as_known  (:pr:`3888`) `Sriharsha Hatwar`_
-  Fix DataFrame.assign for callables (:pr:`3919`) `Tom Augspurger`_
-  Include partitions with no width in repartition (:pr:`3941`) `Matthew Rocklin`_
-  Don't constrict stage/k dtype in dataframe shuffle (:pr:`3942`) `Matthew Rocklin`_

Documentation
^^^^^^^^^^^^^

-  DOC: Add hint on how to render task graphs horizontally (:pr:`3922`) `Uwe Korn`_
-  Add try-now button to main landing page (:pr:`3924`) `Matthew Rocklin`_


.. _v0.19.0 / 2018-08-29:

0.19.0 / 2018-08-29
-------------------

Array
^^^^^

-  Support coordinate in gradient (:pr:`3949`) `Keisuke Fujii`_
-  Fix argtopk split_every bug (:pr:`3810`) `Guido Imperiale`_
-  Ensure result computing dask.array.isnull() always gives a numpy array (:pr:`3825`) `Stephan Hoyer`_
-  Support concatenate for scipy.sparse in dask array (:pr:`3836`) `Matthew Rocklin`_
-  Fix argtopk on 32-bit systems. (:pr:`3823`) `Elliott Sales de Andrade`_
-  Normalize keys in rechunk (:pr:`3820`) `Matthew Rocklin`_
-  Allow shape of dask.array to be a numpy array (:pr:`3844`) `Mark Harfouche`_
-  Fix numpy deprecation warning on tuple indexing (:pr:`3851`) `Tobias de Jong`_
-  Rename ghost module to overlap (:pr:`3830`) `Robert Sare`_
-  Re-add the ghost import to da __init__ (:pr:`3861`) `Jim Crist`_
-  Ensure copy preserves masked arrays (:pr:`3852`) `Tobias de Jong`_

DataFrame
^^^^^^^^^^

-  Added ``dtype`` and ``sparse`` keywords to :func:`dask.dataframe.get_dummies` (:pr:`3792`) `Tom Augspurger`_
-  Added :meth:`dask.dataframe.to_dask_array` for converting a Dask Series or DataFrame to a
   Dask Array, possibly with known chunk sizes (:pr:`3884`) `Tom Augspurger`
-  Changed the behavior for :meth:`dask.array.asarray` for dask dataframe and series inputs. Previously,
   the series was eagerly converted to an in-memory NumPy array before creating a dask array with known
   chunks sizes. This caused unexpectedly high memory usage. Now, no intermediate NumPy array is created,
   and a Dask array with unknown chunk sizes is returned (:pr:`3884`) `Tom Augspurger`
-  DataFrame.iloc (:pr:`3805`) `Tom Augspurger`_
-  When reading multiple paths, expand globs. (:pr:`3828`) `Irina Truong`_
-  Added index column name after resample (:pr:`3833`) `Eric Bonfadini`_
-  Add (lazy) shape property to dataframe and series (:pr:`3212`) `Henrique Ribeiro`_
-  Fix failing hdfs test [test-hdfs] (:pr:`3858`) `Jim Crist`_
-  Fixes for pyarrow 0.10.0 release (:pr:`3860`) `Jim Crist`_
-  Rename to_csv keys for diagnostics (:pr:`3890`) `Matthew Rocklin`_
-  Match pandas warnings for concat sort (:pr:`3897`) `Tom Augspurger`_
-  Include filename in read_csv (:pr:`3908`) `Julia Signell`_

Core
^^^^

-  Better error message on import when missing common dependencies (:pr:`3771`) `Danilo Horta`_
-  Drop Python 3.4 support (:pr:`3840`) `Jim Crist`_
-  Remove expired deprecation warnings (:pr:`3841`) `Jim Crist`_
-  Add DASK_ROOT_CONFIG environment variable (:pr:`3849`) `Joe Hamman`_
-  Don't cull in local scheduler, do cull in delayed (:pr:`3856`) `Jim Crist`_
-  Increase conda download retries (:pr:`3857`) `Jim Crist`_
-  Add python_requires and Trove classifiers (:pr:`3855`) `@hugovk`_
-  Fix collections.abc deprecation warnings in Python 3.7.0 (:pr:`3876`) `Jan Margeta`_
-  Allow dot jpeg to xfail in visualize tests (:pr:`3896`) `Matthew Rocklin`_
-  Add Python 3.7 to travis.yml (:pr:`3894`) `Matthew Rocklin`_
-  Add expand_environment_variables to dask.config (:pr:`3893`) `Joe Hamman`_

Docs
^^^^

-  Fix typo in import statement of diagnostics (:pr:`3826`) `John Mrziglod`_
-  Add link to YARN docs (:pr:`3838`) `Jim Crist`_
-  fix of minor typos in landing page index.html (:pr:`3746`) `Christoph Moehl`_
-  Update delayed-custom.rst (:pr:`3850`) `Anderson Banihirwe`_
-  DOC: clarify delayed docstring (:pr:`3709`) `Scott Sievert`_
-  Add new presentations (:pr:`3880`) `Javad`_
-  Add dask array normalize_chunks to documentation (:pr:`3878`) `Daniel Rothenberg`_
-  Docs: Fix link to snakeviz (:pr:`3900`) `Hans Moritz Günther`_
-  Add missing ` to docstring (:pr:`3915`) `@rtobar`_


.. _v0.18.2 / 2018-07-23:

0.18.2 / 2018-07-23
-------------------

Array
^^^^^

- Reimplemented ``argtopk`` to make it release the GIL (:pr:`3610`) `Guido Imperiale`_
- Don't overlap on non-overlapped dimensions in ``map_overlap`` (:pr:`3653`) `Matthew Rocklin`_
- Fix ``linalg.tsqr`` for dimensions of uncertain length (:pr:`3662`) `Jeremy Chen`_
- Break apart uneven array-of-int slicing to separate chunks (:pr:`3648`) `Matthew Rocklin`_
- Align auto chunks to provided chunks, rather than shape (:pr:`3679`) `Matthew Rocklin`_
- Adds endpoint and retstep support for linspace (:pr:`3675`) `James Bourbeau`_
- Implement ``.blocks`` accessor (:pr:`3689`) `Matthew Rocklin`_
- Add ``block_info`` keyword to ``map_blocks`` functions (:pr:`3686`) `Matthew Rocklin`_
- Slice by dask array of ints (:pr:`3407`) `Guido Imperiale`_
- Support ``dtype`` in ``arange`` (:pr:`3722`) `Guido Imperiale`_
- Fix ``argtopk`` with uneven chunks (:pr:`3720`) `Guido Imperiale`_
- Raise error when ``replace=False`` in ``da.choice`` (:pr:`3765`) `James Bourbeau`_
- Update chunks in ``Array.__setitem__`` (:pr:`3767`) `Itamar Turner-Trauring`_
- Add a ``chunksize`` convenience property (:pr:`3777`) `Jacob Tomlinson`_
- Fix and simplify array slicing behavior when ``step < 0`` (:pr:`3702`) `Ziyao Wei`_
- Ensure ``to_zarr`` with ``return_stored`` ``True`` returns a Dask Array (:pr:`3786`) `John A Kirkham`_

Bag
^^^

- Add ``last_endline`` optional parameter in ``to_textfiles`` (:pr:`3745`) `George Sakkis`_

Dataframe
^^^^^^^^^

- Add aggregate function for rolling objects (:pr:`3772`) `Gerome Pistre`_
- Properly tokenize cumulative groupby aggregations (:pr:`3799`) `Cloves Almeida`_

Delayed
^^^^^^^

- Add the ``@`` operator to the delayed objects (:pr:`3691`) `Mark Harfouche`_
- Add delayed best practices to documentation (:pr:`3737`) `Matthew Rocklin`_
- Fix ``@delayed`` decorator for methods and add tests (:pr:`3757`) `Ziyao Wei`_

Core
^^^^

- Fix extra progressbar (:pr:`3669`) `Mike Neish`_
- Allow tasks back onto ordering stack if they have one dependency (:pr:`3652`) `Matthew Rocklin`_
- Prefer end-tasks with low numbers of dependencies when ordering (:pr:`3588`) `Tom Augspurger`_
- Add ``assert_eq`` to top-level modules (:pr:`3726`) `Matthew Rocklin`_
- Test that dask collections can hold ``scipy.sparse`` arrays (:pr:`3738`) `Matthew Rocklin`_
- Fix setup of lz4 decompression functions (:pr:`3782`) `Elliott Sales de Andrade`_
- Add datasets module (:pr:`3780`) `Matthew Rocklin`_


.. _v0.18.1 / 2018-06-22:

0.18.1 / 2018-06-22
-------------------

Array
^^^^^

- ``from_array`` now supports scalar types and nested lists/tuples in input,
  just like all numpy functions do; it also produces a simpler graph when the
  input is a plain ndarray (:pr:`3568`) `Guido Imperiale`_
- Fix slicing of big arrays due to cumsum dtype bug (:pr:`3620`) `Marco Rossi`_
- Add Dask Array implementation of pad (:pr:`3578`) `John A Kirkham`_
- Fix array random API examples (:pr:`3625`) `James Bourbeau`_
- Add average function to dask array (:pr:`3640`) `James Bourbeau`_
- Tokenize ghost_internal with axes (:pr:`3643`)  `Matthew Rocklin`_
- Add outer for Dask Arrays (:pr:`3658`) `John A Kirkham`_

DataFrame
^^^^^^^^^

- Add Index.to_series method (:pr:`3613`) `Henrique Ribeiro`_
- Fix missing partition columns in pyarrow-parquet (:pr:`3636`) `Martin Durant`_

Core
^^^^

- Minor tweaks to CI (:pr:`3629`) `Guido Imperiale`_
- Add back dask.utils.effective_get (:pr:`3642`) `Matthew Rocklin`_
- DASK_CONFIG dictates config write location (:pr:`3621`) `Jim Crist`_
- Replace 'collections' key in unpack_collections with unique key (:pr:`3632`) `Yu Feng`_
- Avoid deepcopy in dask.config.set (:pr:`3649`) `Matthew Rocklin`_


.. _v0.18.0 / 2018-06-14:

0.18.0 / 2018-06-14
-------------------

Array
^^^^^

- Add to/from_zarr for Zarr-format datasets and arrays (:pr:`3460`) `Martin Durant`_
- Experimental addition of generalized ufunc support, ``apply_gufunc``, ``gufunc``, and
  ``as_gufunc`` (:pr:`3109`) (:pr:`3526`) (:pr:`3539`) `Markus Gonser`_
- Avoid unnecessary rechunking tasks (:pr:`3529`) `Matthew Rocklin`_
- Compute dtypes at runtime for fft (:pr:`3511`) `Matthew Rocklin`_
- Generate UUIDs for all da.store operations (:pr:`3540`) `Martin Durant`_
- Correct internal dimension of Dask's SVD (:pr:`3517`) `John A Kirkham`_
- BUG: do not raise IndexError for identity slice in array.vindex (:pr:`3559`) `Scott Sievert`_
- Adds `isneginf` and `isposinf` (:pr:`3581`) `John A Kirkham`_
- Drop Dask Array's `learn` module (:pr:`3580`) `John A Kirkham`_
- added sfqr (short-and-fat) as a counterpart to tsqr… (:pr:`3575`) `Jeremy Chen`_
- Allow 0-width chunks in dask.array.rechunk (:pr:`3591`) `Marc Pfister`_
- Document Dask Array's `nan_to_num` in public API (:pr:`3599`) `John A Kirkham`_
- Show block example (:pr:`3601`) `John A Kirkham`_
- Replace token= keyword with name= in map_blocks (:pr:`3597`) `Matthew Rocklin`_
- Disable locking in to_zarr (needed for using to_zarr in a distributed context) (:pr:`3607`) `John A Kirkham`_
- Support Zarr Arrays in `to_zarr`/`from_zarr` (:pr:`3561`) `John A Kirkham`_
- Added recursion to array/linalg/tsqr to better manage the single core bottleneck (:pr:`3586`) `Jeremy Chan`_
  (:pr:`3396`) `Guido Imperiale`_

Dataframe
^^^^^^^^^

- Add to/read_json (:pr:`3494`) `Martin Durant`_
- Adds ``index`` to unsupported arguments for ``DataFrame.rename`` method (:pr:`3522`) `James Bourbeau`_
- Adds support to subset Dask DataFrame columns using ``numpy.ndarray``, ``pandas.Series``, and
  ``pandas.Index`` objects (:pr:`3536`) `James Bourbeau`_
- Raise error if meta columns do not match dataframe (:pr:`3485`) `Christopher Ren`_
- Add index to unsupprted argument for DataFrame.rename (:pr:`3522`) `James Bourbeau`_
- Adds support for subsetting DataFrames with pandas Index/Series and numpy ndarrays (:pr:`3536`) `James Bourbeau`_
- Dataframe sample method docstring fix (:pr:`3566`) `James Bourbeau`_
- fixes dd.read_json to infer file compression (:pr:`3594`) `Matt Lee`_
- Adds n to sample method (:pr:`3606`) `James Bourbeau`_
- Add fastparquet ParquetFile object support (:pr:`3573`) `@andrethrill`_

Bag
^^^

- Rename method= keyword to shuffle= in bag.groupby (:pr:`3470`) `Matthew Rocklin`_

Core
^^^^

- Replace get= keyword with scheduler= keyword (:pr:`3448`) `Matthew Rocklin`_
- Add centralized dask.config module to handle configuration for all Dask
  subprojects (:pr:`3432`) (:pr:`3513`) (:pr:`3520`) `Matthew Rocklin`_
- Add `dask-ssh` CLI Options and Description. (:pr:`3476`) `@beomi`_
- Read whole files fix regardless of header for HTTP (:pr:`3496`) `Martin Durant`_
- Adds synchronous scheduler syntax to debugging docs (:pr:`3509`) `James Bourbeau`_
- Replace dask.set_options with dask.config.set (:pr:`3502`) `Matthew Rocklin`_
- Update sphinx readthedocs-theme (:pr:`3516`) `Matthew Rocklin`_
- Introduce "auto" value for normalize_chunks (:pr:`3507`) `Matthew Rocklin`_
- Fix check in configuration with env=None (:pr:`3562`) `Simon Perkins`_
- Update sizeof definitions (:pr:`3582`) `Matthew Rocklin`_
- Remove --verbose flag from travis-ci (:pr:`3477`) `Matthew Rocklin`_
- Remove "da.random" from random array keys (:pr:`3604`) `Matthew Rocklin`_


.. _v0.17.5 / 2018-05-16:

0.17.5 / 2018-05-16
-------------------

Array
^^^^^

- Fix ``rechunk`` with chunksize of -1 in a dict (:pr:`3469`) `Stephan Hoyer`_
- ``einsum`` now accepts the ``split_every`` parameter (:pr:`3471`) `Guido Imperiale`_
- Improved slicing performance (:pr:`3479`) `Yu Feng`_

DataFrame
^^^^^^^^^

- Compatibility with pandas 0.23.0 (:pr:`3499`) `Tom Augspurger`_


.. _v0.17.4 / 2018-05-03:

0.17.4 / 2018-05-03
-------------------

Dataframe
^^^^^^^^^

- Add support for indexing Dask DataFrames with string subclasses (:pr:`3461`) `James Bourbeau`_
- Allow using both sorted_index and chunksize in read_hdf (:pr:`3463`) `Pierre Bartet`_
- Pass filesystem to arrow piece reader (:pr:`3466`) `Martin Durant`_
- Switches to using dask.compat string_types (:pr:`3462`) `James Bourbeau`_


.. _v0.17.3 / 2018-05-02:

0.17.3 / 2018-05-02
-------------------

Array
^^^^^

- Add ``einsum`` for Dask Arrays (:pr:`3412`) `Simon Perkins`_
- Add ``piecewise`` for Dask Arrays (:pr:`3350`) `John A Kirkham`_
- Fix handling of ``nan`` in ``broadcast_shapes`` (:pr:`3356`) `John A Kirkham`_
- Add ``isin`` for dask arrays (:pr:`3363`). `Stephan Hoyer`_
- Overhauled ``topk`` for Dask Arrays: faster algorithm, particularly for large k's; added support
  for multiple axes, recursive aggregation, and an option to pick the bottom k elements instead.
  (:pr:`3395`) `Guido Imperiale`_
- The ``topk`` API has changed from topk(k, array) to the more conventional topk(array, k).
  The legacy API still works but is now deprecated. (:pr:`2965`) `Guido Imperiale`_
- New function ``argtopk`` for Dask Arrays (:pr:`3396`) `Guido Imperiale`_
- Fix handling partial depth and boundary in ``map_overlap`` (:pr:`3445`) `John A Kirkham`_
- Add ``gradient`` for Dask Arrays (:pr:`3434`) `John A Kirkham`_

DataFrame
^^^^^^^^^

- Allow `t` as shorthand for `table` in `to_hdf` for pandas compatibility (:pr:`3330`) `Jörg Dietrich`_
- Added top level `isna` method for Dask DataFrames (:pr:`3294`) `Christopher Ren`_
- Fix selection on partition column on ``read_parquet`` for ``engine="pyarrow"`` (:pr:`3207`) `Uwe Korn`_
- Added DataFrame.squeeze method (:pr:`3366`) `Christopher Ren`_
- Added `infer_divisions` option to ``read_parquet`` to specify whether read engines should compute divisions (:pr:`3387`) `Jon Mease`_
- Added support for inferring division for ``engine="pyarrow"`` (:pr:`3387`) `Jon Mease`_
- Provide more informative error message for meta= errors (:pr:`3343`) `Matthew Rocklin`_
- add orc reader (:pr:`3284`) `Martin Durant`_
- Default compression for parquet now always Snappy, in line with pandas (:pr:`3373`) `Martin Durant`_
- Fixed bug in Dask DataFrame and Series comparisons with NumPy scalars (:pr:`3436`) `James Bourbeau`_
- Remove outdated requirement from repartition docstring (:pr:`3440`) `Jörg Dietrich`_
- Fixed bug in aggregation when only a Series is selected (:pr:`3446`) `Jörg Dietrich`_
- Add default values to make_timeseries (:pr:`3421`) `Matthew Rocklin`_

Core
^^^^

- Support traversing collections in persist, visualize, and optimize (:pr:`3410`) `Jim Crist`_
- Add schedule= keyword to compute and persist.  This replaces common use of the get= keyword (:pr:`3448`) `Matthew Rocklin`_


.. _v0.17.2 / 2018-03-21:

0.17.2 / 2018-03-21
-------------------

Array
^^^^^

- Add ``broadcast_arrays`` for Dask Arrays (:pr:`3217`) `John A Kirkham`_
- Add ``bitwise_*`` ufuncs (:pr:`3219`) `John A Kirkham`_
- Add optional ``axis`` argument to ``squeeze`` (:pr:`3261`) `John A Kirkham`_
- Validate inputs to atop (:pr:`3307`) `Matthew Rocklin`_
- Avoid calls to astype in concatenate if all parts have the same dtype (:pr:`3301`) `Martin Durant`_

DataFrame
^^^^^^^^^

- Fixed bug in shuffle due to aggressive truncation (:pr:`3201`) `Matthew Rocklin`_
- Support specifying categorical columns on ``read_parquet`` with ``categories=[…]`` for ``engine="pyarrow"`` (:pr:`3177`) `Uwe Korn`_
- Add ``dd.tseries.Resampler.agg`` (:pr:`3202`) `Richard Postelnik`_
- Support operations that mix dataframes and arrays (:pr:`3230`) `Matthew Rocklin`_
- Support extra Scalar and Delayed args in ``dd.groupby._Groupby.apply`` (:pr:`3256`) `Gabriele Lanaro`_

Bag
^^^

- Support joining against single-partitioned bags and delayed objects (:pr:`3254`) `Matthew Rocklin`_

Core
^^^^

- Fixed bug when using unexpected but hashable types for keys (:pr:`3238`) `Daniel Collins`_
- Fix bug in task ordering so that we break ties consistently with the key name (:pr:`3271`) `Matthew Rocklin`_
- Avoid sorting tasks in order when the number of tasks is very large (:pr:`3298`) `Matthew Rocklin`_


.. _v0.17.1 / 2018-02-22:

0.17.1 / 2018-02-22
-------------------

Array
^^^^^

- Corrected dimension chunking in indices (:issue:`3166`, :pr:`3167`) `Simon Perkins`_
- Inline ``store_chunk`` calls for ``store``'s ``return_stored`` option (:pr:`3153`) `John A Kirkham`_
- Compatibility with struct dtypes for NumPy 1.14.1 release (:pr:`3187`) `Matthew Rocklin`_

DataFrame
^^^^^^^^^

- Bugfix to allow column assignment of pandas datetimes(:pr:`3164`) `Max Epstein`_

Core
^^^^

- New file-system for HTTP(S), allowing direct loading from specific URLs (:pr:`3160`) `Martin Durant`_
- Fix bug when tokenizing partials with no keywords (:pr:`3191`) `Matthew Rocklin`_
- Use more recent LZ4 API (:pr:`3157`) `Thrasibule`_
- Introduce output stream parameter for progress bar (:pr:`3185`) `Dieter Weber`_


.. _v0.17.0 / 2018-02-09:

0.17.0 / 2018-02-09
-------------------

Array
^^^^^

- Added a support object-type arrays for nansum, nanmin, and nanmax (:issue:`3133`) `Keisuke Fujii`_
- Update error handling when len is called with empty chunks (:issue:`3058`) `Xander Johnson`_
- Fixes a metadata bug with ``store``'s ``return_stored`` option (:pr:`3064`) `John A Kirkham`_
- Fix a bug in ``optimization.fuse_slice`` to properly handle when first input is ``None`` (:pr:`3076`) `James Bourbeau`_
- Support arrays with unknown chunk sizes in percentile (:pr:`3107`) `Matthew Rocklin`_
- Tokenize scipy.sparse arrays and np.matrix (:pr:`3060`) `Roman Yurchak`_

DataFrame
^^^^^^^^^
- Support month timedeltas in repartition(freq=...) (:pr:`3110`) `Matthew Rocklin`_
- Avoid mutation in dataframe groupby tests (:pr:`3118`) `Matthew Rocklin`_
- ``read_csv``, ``read_table``, and ``read_parquet`` accept iterables of paths
  (:pr:`3124`) `Jim Crist`_
- Deprecates the ``dd.to_delayed`` *function* in favor of the existing method
  (:pr:`3126`) `Jim Crist`_
- Return dask.arrays from df.map_partitions calls when the UDF returns a numpy array (:pr:`3147`) `Matthew Rocklin`_
- Change handling of ``columns`` and ``index`` in ``dd.read_parquet`` to be more
  consistent, especially in handling of multi-indices (:pr:`3149`) `Jim Crist`_
- fastparquet append=True allowed to create new dataset (:pr:`3097`) `Martin Durant`_
- dtype rationalization for sql queries (:pr:`3100`) `Martin Durant`_

Bag
^^^

- Document ``bag.map_paritions`` function may receive either a list or generator. (:pr:`3150`) `Nir`_

Core
^^^^

- Change default task ordering to prefer nodes with few dependents and then
  many downstream dependencies (:pr:`3056`) `Matthew Rocklin`_
- Add color= option to visualize to color by task order (:pr:`3057`) (:pr:`3122`) `Matthew Rocklin`_
- Deprecate ``dask.bytes.open_text_files`` (:pr:`3077`) `Jim Crist`_
- Remove short-circuit hdfs reads handling due to maintenance costs. May be
  re-added in a more robust manner later (:pr:`3079`) `Jim Crist`_
- Add ``dask.base.optimize`` for optimizing multiple collections without
  computing. (:pr:`3071`) `Jim Crist`_
- Rename ``dask.optimize`` module to ``dask.optimization`` (:pr:`3071`) `Jim Crist`_
- Change task ordering to do a full traversal (:pr:`3066`) `Matthew Rocklin`_
- Adds an ``optimize_graph`` keyword to all ``to_delayed`` methods to allow
  controlling whether optimizations occur on conversion. (:pr:`3126`) `Jim Crist`_
- Support using ``pyarrow`` for hdfs integration (:pr:`3123`) `Jim Crist`_
- Move HDFS integration and tests into dask repo (:pr:`3083`) `Jim Crist`_
- Remove write_bytes (:pr:`3116`) `Jim Crist`_


.. _v0.16.1 / 2018-01-09:

0.16.1 / 2018-01-09
-------------------

Array
^^^^^

- Fix handling of scalar percentile values in ``percentile`` (:pr:`3021`) `James Bourbeau`_
- Prevent ``bool()`` coercion from calling compute (:pr:`2958`) `Albert DeFusco`_
- Add ``matmul`` (:pr:`2904`) `John A Kirkham`_
- Support N-D arrays with ``matmul`` (:pr:`2909`) `John A Kirkham`_
- Add ``vdot`` (:pr:`2910`) `John A Kirkham`_
- Explicit ``chunks`` argument for ``broadcast_to`` (:pr:`2943`) `Stephan Hoyer`_
- Add ``meshgrid`` (:pr:`2938`) `John A Kirkham`_ and (:pr:`3001`) `Markus Gonser`_
- Preserve singleton chunks in ``fftshift``/``ifftshift`` (:pr:`2733`) `John A Kirkham`_
- Fix handling of negative indexes in ``vindex`` and raise errors for out of bounds indexes (:pr:`2967`) `Stephan Hoyer`_
- Add ``flip``, ``flipud``, ``fliplr`` (:pr:`2954`) `John A Kirkham`_
- Add ``float_power`` ufunc (:pr:`2962`) (:pr:`2969`) `John A Kirkham`_
- Compatibility for changes to structured arrays in the upcoming NumPy 1.14 release (:pr:`2964`) `Tom Augspurger`_
- Add ``block`` (:pr:`2650`) `John A Kirkham`_
- Add ``frompyfunc`` (:pr:`3030`) `Jim Crist`_
- Add the ``return_stored`` option to ``store`` for chaining stored results (:pr:`2980`) `John A Kirkham`_

DataFrame
^^^^^^^^^

- Fixed naming bug in cumulative aggregations (:issue:`3037`) `Martijn Arts`_
- Fixed ``dd.read_csv`` when ``names`` is given but ``header`` is not set to ``None`` (:issue:`2976`) `Martijn Arts`_
- Fixed ``dd.read_csv`` so that passing instances of ``CategoricalDtype`` in ``dtype`` will result in known categoricals (:pr:`2997`) `Tom Augspurger`_
- Prevent ``bool()`` coercion from calling compute (:pr:`2958`) `Albert DeFusco`_
- ``DataFrame.read_sql()`` (:pr:`2928`) to an empty database tables returns an empty dask dataframe `Apostolos Vlachopoulos`_
- Compatibility for reading Parquet files written by PyArrow 0.8.0 (:pr:`2973`) `Tom Augspurger`_
- Correctly handle the column name (`df.columns.name`) when reading in ``dd.read_parquet`` (:pr:`2973`) `Tom Augspurger`_
- Fixed ``dd.concat`` losing the index dtype when the data contained a categorical (:issue:`2932`) `Tom Augspurger`_
- Add ``dd.Series.rename`` (:pr:`3027`) `Jim Crist`_
- ``DataFrame.merge()`` now supports merging on a combination of columns and the index (:pr:`2960`) `Jon Mease`_
- Removed the deprecated ``dd.rolling*`` methods, in preparation for their removal in the next pandas release (:pr:`2995`) `Tom Augspurger`_
- Fix metadata inference bug in which single-partition series were mistakenly special cased (:pr:`3035`) `Jim Crist`_
- Add support for ``Series.str.cat`` (:pr:`3028`) `Jim Crist`_

Core
^^^^

- Improve 32-bit compatibility (:pr:`2937`) `Matthew Rocklin`_
- Change task prioritization to avoid upwards branching (:pr:`3017`) `Matthew Rocklin`_


.. _v0.16.0 / 2017-11-17:

0.16.0 / 2017-11-17
-------------------

This is a major release.  It includes breaking changes, new protocols, and a
large number of bug fixes.

Array
^^^^^

- Add ``atleast_1d``, ``atleast_2d``, and ``atleast_3d`` (:pr:`2760`) (:pr:`2765`) `John A Kirkham`_
- Add ``allclose`` (:pr:`2771`) by `John A Kirkham`_
- Remove ``random.different_seeds`` from Dask Array API docs (:pr:`2772`) `John A Kirkham`_
- Deprecate ``vnorm`` in favor of ``dask.array.linalg.norm`` (:pr:`2773`) `John A Kirkham`_
- Reimplement ``unique`` to be lazy (:pr:`2775`) `John A Kirkham`_
- Support broadcasting of Dask Arrays with 0-length dimensions (:pr:`2784`) `John A Kirkham`_
- Add ``asarray`` and ``asanyarray`` to Dask Array API docs (:pr:`2787`) `James Bourbeau`_
- Support ``unique``'s ``return_*`` arguments (:pr:`2779`) `John A Kirkham`_
- Simplify ``_unique_internal`` (:pr:`2850`) (:pr:`2855`) `John A Kirkham`_
- Avoid removing some getter calls in array optimizations (:pr:`2826`) `Jim Crist`_

DataFrame
^^^^^^^^^

- Support ``pyarrow`` in ``dd.to_parquet`` (:pr:`2868`) `Jim Crist`_
- Fixed ``DataFrame.quantile`` and ``Series.quantile`` returning ``nan`` when missing values are present (:pr:`2791`) `Tom Augspurger`_
- Fixed ``DataFrame.quantile`` losing the result ``.name`` when ``q`` is a scalar (:pr:`2791`) `Tom Augspurger`_
- Fixed ``dd.concat`` return a ``dask.Dataframe`` when concatenating a single series along the columns, matching pandas' behavior (:pr:`2800`) `James Munroe`_
- Fixed default inplace parameter for ``DataFrame.eval`` to match the pandas defualt for pandas >= 0.21.0 (:pr:`2838`) `Tom Augspurger`_
- Fix exception when calling ``DataFrame.set_index`` on text column where one of the partitions was empty (:pr:`2831`) `Jesse Vogt`_
- Do not raise exception when calling ``DataFrame.set_index`` on empty dataframe (:pr:`2827`) `Jesse Vogt`_
- Fixed bug in ``Dataframe.fillna`` when filling with a ``Series`` value (:pr:`2810`) `Tom Augspurger`_
- Deprecate old argument ordering in ``dd.to_parquet`` to better match convention of putting the dataframe first (:pr:`2867`) `Jim Crist`_
- df.astype(categorical_dtype -> known categoricals (:pr:`2835`) `Jim Crist`_
- Test against Pandas release candidate (:pr:`2814`) `Tom Augspurger`_
- Add more tests for read_parquet(engine='pyarrow') (:pr:`2822`) `Uwe Korn`_
- Remove unnecessary map_partitions in aggregate (:pr:`2712`) `Christopher Prohm`_
- Fix bug calling sample on empty partitions (:pr:`2818`) `@xwang777`_
- Error nicely when parsing dates in read_csv (:pr:`2863`) `Jim Crist`_
- Cleanup handling of passing filesystem objects to PyArrow readers (:pr:`2527`) `@fjetter`_
- Support repartitioning even if there are no divisions (:pr:`2873`) `@Ced4`_
- Support reading/writing to hdfs using ``pyarrow`` in ``dd.to_parquet`` (:pr:`2894`, :pr:`2881`) `Jim Crist`_

Core
^^^^

-  Allow tuples as sharedict keys (:pr:`2763`) `Matthew Rocklin`_
-  Calling compute within a dask.distributed task defaults to distributed scheduler (:pr:`2762`) `Matthew Rocklin`_
-  Auto-import gcsfs when gcs:// protocol is used (:pr:`2776`) `Matthew Rocklin`_
-  Fully remove dask.async module, use dask.local instead (:pr:`2828`) `Thomas Caswell`_
-  Compatibility with bokeh 0.12.10 (:pr:`2844`) `Tom Augspurger`_
-  Reduce test memory usage (:pr:`2782`) `Jim Crist`_
-  Add Dask collection interface (:pr:`2748`) `Jim Crist`_
-  Update Dask collection interface during XArray integration (:pr:`2847`) `Matthew Rocklin`_
-  Close resource profiler process on __exit__ (:pr:`2871`) `Jim Crist`_
-  Fix S3 tests (:pr:`2875`) `Jim Crist`_
-  Fix port for bokeh dashboard in docs (:pr:`2889`) `Ian Hopkinson`_
-  Wrap Dask filesystems for PyArrow compatibility (:pr:`2881`) `Jim Crist`_


.. _v0.15.4 / 2017-10-06:

0.15.4 / 2017-10-06
-------------------

Array
^^^^^

-  ``da.random.choice`` now works with array arguments (:pr:`2781`)
-  Support indexing in arrays with np.int (fixes regression) (:pr:`2719`)
-  Handle zero dimension with rechunking (:pr:`2747`)
-  Support -1 as an alias for "size of the dimension" in ``chunks`` (:pr:`2749`)
-  Call mkdir in array.to_npy_stack (:pr:`2709`)

DataFrame
^^^^^^^^^

-  Added the `.str` accessor to Categoricals with string categories (:pr:`2743`)
-  Support int96 (spark) datetimes in parquet writer (:pr:`2711`)
-  Pass on file scheme to fastparquet (:pr:`2714`)
-  Support Pandas 0.21 (:pr:`2737`)

Bag
^^^

- Add tree reduction support for foldby (:pr:`2710`)


Core
^^^^

-  Drop s3fs from ``pip install dask[complete]`` (:pr:`2750`)


.. _v0.15.3 / 2017-09-24:

0.15.3 / 2017-09-24
-------------------

Array
^^^^^

-  Add masked arrays (:pr:`2301`)
-  Add ``*_like array creation functions`` (:pr:`2640`)
-  Indexing with unsigned integer array (:pr:`2647`)
-  Improved slicing with boolean arrays of different dimensions (:pr:`2658`)
-  Support literals in ``top`` and ``atop`` (:pr:`2661`)
-  Optional axis argument in cumulative functions (:pr:`2664`)
-  Improve tests on scalars with ``assert_eq`` (:pr:`2681`)
-  Fix norm keepdims (:pr:`2683`)
-  Add ``ptp`` (:pr:`2691`)
-  Add apply_along_axis (:pr:`2690`) and apply_over_axes (:pr:`2702`)

DataFrame
^^^^^^^^^

-  Added ``Series.str[index]`` (:pr:`2634`)
-  Allow the groupby by param to handle columns and index levels (:pr:`2636`)
-  ``DataFrame.to_csv`` and ``Bag.to_textfiles`` now return the filenames to
    which they have written (:pr:`2655`)
-  Fix combination of ``partition_on`` and ``append`` in ``to_parquet``
   (:pr:`2645`)
-  Fix for parquet file schemes (:pr:`2667`)
-  Repartition works with mixed categoricals (:pr:`2676`)

Core
^^^^

-  ``python setup.py test`` now runs tests (:pr:`2641`)
-  Added new cheatsheet (:pr:`2649`)
-  Remove resize tool in Bokeh plots (:pr:`2688`)


.. _v0.15.2 / 2017-08-25:

0.15.2 / 2017-08-25
-------------------

Array
^^^^^

-  Remove spurious keys from map_overlap graph (:pr:`2520`)
-  where works with non-bool condition and scalar values (:pr:`2543`) (:pr:`2549`)
-  Improve compress (:pr:`2541`) (:pr:`2545`) (:pr:`2555`)
-  Add argwhere, _nonzero, and where(cond) (:pr:`2539`)
-  Generalize vindex in dask.array to handle multi-dimensional indices (:pr:`2573`)
-  Add choose method (:pr:`2584`)
-  Split code into reorganized files (:pr:`2595`)
-  Add linalg.norm (:pr:`2597`)
-  Add diff, ediff1d (:pr:`2607`), (:pr:`2609`)
-  Improve dtype inference and reflection (:pr:`2571`)

Bag
^^^

-   Remove deprecated Bag behaviors (:pr:`2525`)

DataFrame
^^^^^^^^^

-  Support callables in assign (:pr:`2513`)
-  better error messages for read_csv (:pr:`2522`)
-  Add dd.to_timedelta (:pr:`2523`)
-  Verify metadata in from_delayed (:pr:`2534`) (:pr:`2591`)
-  Add DataFrame.isin (:pr:`2558`)
-  Read_hdf supports iterables of files (:pr:`2547`)

Core
^^^^

-  Remove bare ``except:`` blocks everywhere (:pr:`2590`)

.. _v0.15.1 / 2017-07-08:

0.15.1 / 2017-07-08
-------------------

-  Add storage_options to to_textfiles and to_csv (:pr:`2466`)
-  Rechunk and simplify rfftfreq (:pr:`2473`), (:pr:`2475`)
-  Better support ndarray subclasses (:pr:`2486`)
-  Import star in dask.distributed (:pr:`2503`)
-  Threadsafe cache handling with tokenization (:pr:`2511`)


.. _v0.15.0 / 2017-06-09:

0.15.0 / 2017-06-09
-------------------

Array
^^^^^

-  Add dask.array.stats submodule (:pr:`2269`)
-  Support ``ufunc.outer`` (:pr:`2345`)
-  Optimize fancy indexing by reducing graph overhead (:pr:`2333`) (:pr:`2394`)
-  Faster array tokenization using alternative hashes (:pr:`2377`)
-  Added the matmul ``@`` operator (:pr:`2349`)
-  Improved coverage of the ``numpy.fft`` module (:pr:`2320`) (:pr:`2322`) (:pr:`2327`) (:pr:`2323`)
-  Support NumPy's ``__array_ufunc__`` protocol (:pr:`2438`)

Bag
^^^

-  Fix bug where reductions on bags with no partitions would fail (:pr:`2324`)
-  Add broadcasting and variadic ``db.map`` top-level function.  Also remove
   auto-expansion of tuples as map arguments (:pr:`2339`)
-  Rename ``Bag.concat`` to ``Bag.flatten`` (:pr:`2402`)

DataFrame
^^^^^^^^^

-  Parquet improvements (:pr:`2277`) (:pr:`2422`)

Core
^^^^

-  Move dask.async module to dask.local (:pr:`2318`)
-  Support callbacks with nested scheduler calls (:pr:`2397`)
-  Support pathlib.Path objects as uris  (:pr:`2310`)


.. _v0.14.3 / 2017-05-05:

0.14.3 / 2017-05-05
-------------------

DataFrame
^^^^^^^^^

-  Pandas 0.20.0 support

.. _v0.14.2 / 2017-05-03:

0.14.2 / 2017-05-03
-------------------

Array
^^^^^

-  Add da.indices (:pr:`2268`), da.tile (:pr:`2153`), da.roll (:pr:`2135`)
-  Simultaneously support drop_axis and new_axis in da.map_blocks (:pr:`2264`)
-  Rechunk and concatenate work with unknown chunksizes (:pr:`2235`) and (:pr:`2251`)
-  Support non-numpy container arrays, notably sparse arrays (:pr:`2234`)
-  Tensordot contracts over multiple axes (:pr:`2186`)
-  Allow delayed targets in da.store (:pr:`2181`)
-  Support interactions against lists and tuples (:pr:`2148`)
-  Constructor plugins for debugging (:pr:`2142`)
-  Multi-dimensional FFTs (single chunk) (:pr:`2116`)

Bag
^^^

-  to_dataframe enforces consistent types (:pr:`2199`)

DataFrame
^^^^^^^^^

-  Set_index always fully sorts the index (:pr:`2290`)
-  Support compatibility with pandas 0.20.0 (:pr:`2249`), (:pr:`2248`), and (:pr:`2246`)
-  Support Arrow Parquet reader (:pr:`2223`)
-  Time-based rolling windows (:pr:`2198`)
-  Repartition can now create more partitions, not just less (:pr:`2168`)

Core
^^^^

-  Always use absolute paths when on POSIX file system (:pr:`2263`)
-  Support user provided graph optimizations (:pr:`2219`)
-  Refactor path handling (:pr:`2207`)
-  Improve fusion performance (:pr:`2129`), (:pr:`2131`), and (:pr:`2112`)


.. _v0.14.1 / 2017-03-22:

0.14.1 / 2017-03-22
-------------------

Array
^^^^^

-  Micro-optimize optimizations (:pr:`2058`)
-  Change slicing optimizations to avoid fusing raw numpy arrays (:pr:`2075`) (:pr:`2080`)
-  Dask.array operations now work on numpy arrays (:pr:`2079`)
-  Reshape now works in a much broader set of cases (:pr:`2089`)
-  Support deepcopy python protocol (:pr:`2090`)
-  Allow user-provided FFT implementations in ``da.fft`` (:pr:`2093`)

DataFrame
^^^^^^^^^

-  Fix to_parquet with empty partitions (:pr:`2020`)
-  Optional ``npartitions='auto'`` mode in ``set_index`` (:pr:`2025`)
-  Optimize shuffle performance (:pr:`2032`)
-  Support efficient repartitioning along time windows like ``repartition(freq='12h')`` (:pr:`2059`)
-  Improve speed of categorize (:pr:`2010`)
-  Support single-row dataframe arithmetic (:pr:`2085`)
-  Automatically avoid shuffle when setting index with a sorted column (:pr:`2091`)
-  Improve handling of integer-na handling in read_csv (:pr:`2098`)

Delayed
^^^^^^^

-  Repeated attribute access on delayed objects uses the same key (:pr:`2084`)

Core
^^^^

-   Improve naming of nodes in dot visuals to avoid generic ``apply``
    (:pr:`2070`)
-   Ensure that worker processes have different random seeds (:pr:`2094`)


.. _v0.14.0 / 2017-02-24:

0.14.0 / 2017-02-24
-------------------

Array
^^^^^

- Fix corner cases with zero shape and misaligned values in ``arange``
  (:pr:`1902`), (:pr:`1904`), (:pr:`1935`), (:pr:`1955`), (:pr:`1956`)
- Improve concatenation efficiency (:pr:`1923`)
- Avoid hashing in ``from_array`` if name is provided (:pr:`1972`)

Bag
^^^

- Repartition can now increase number of partitions (:pr:`1934`)
- Fix bugs in some reductions with empty partitions (:pr:`1939`), (:pr:`1950`),
  (:pr:`1953`)


DataFrame
^^^^^^^^^

- Support non-uniform categoricals (:pr:`1877`), (:pr:`1930`)
- Groupby cumulative reductions (:pr:`1909`)
- DataFrame.loc indexing now supports lists (:pr:`1913`)
- Improve multi-level groupbys (:pr:`1914`)
- Improved HTML and string repr for DataFrames (:pr:`1637`)
- Parquet append (:pr:`1940`)
- Add ``dd.demo.daily_stock`` function for teaching (:pr:`1992`)

Delayed
^^^^^^^

- Add ``traverse=`` keyword to delayed to optionally avoid traversing nested
  data structures (:pr:`1899`)
- Support Futures in from_delayed functions (:pr:`1961`)
- Improve serialization of decorated delayed functions (:pr:`1969`)

Core
^^^^

- Improve windows path parsing in corner cases (:pr:`1910`)
- Rename tasks when fusing (:pr:`1919`)
- Add top level ``persist`` function (:pr:`1927`)
- Propagate ``errors=`` keyword in byte handling (:pr:`1954`)
- Dask.compute traverses Python collections (:pr:`1975`)
- Structural sharing between graphs in dask.array and dask.delayed (:pr:`1985`)


.. _v0.13.0 / 2017-01-02:

0.13.0 / 2017-01-02
-------------------

Array
^^^^^

- Mandatory dtypes on dask.array.  All operations maintain dtype information
  and UDF functions like map_blocks now require a dtype= keyword if it can not
  be inferred.  (:pr:`1755`)
- Support arrays without known shapes, such as arises when slicing arrays with
  arrays or converting dataframes to arrays (:pr:`1838`)
- Support mutation by setting one array with another (:pr:`1840`)
- Tree reductions for covariance and correlations.  (:pr:`1758`)
- Add SerializableLock for better use with distributed scheduling (:pr:`1766`)
- Improved atop support (:pr:`1800`)
- Rechunk optimization (:pr:`1737`), (:pr:`1827`)

Bag
^^^

- Avoid wrong results when recomputing the same groupby twice (:pr:`1867`)

DataFrame
^^^^^^^^^

- Add ``map_overlap`` for custom rolling operations (:pr:`1769`)
- Add ``shift`` (:pr:`1773`)
- Add Parquet support (:pr:`1782`) (:pr:`1792`) (:pr:`1810`), (:pr:`1843`),
  (:pr:`1859`), (:pr:`1863`)
- Add missing methods combine, abs, autocorr, sem, nsmallest, first, last,
  prod, (:pr:`1787`)
- Approximate nunique (:pr:`1807`), (:pr:`1824`)
- Reductions with multiple output partitions (for operations like
  drop_duplicates) (:pr:`1808`), (:pr:`1823`) (:pr:`1828`)
- Add delitem and copy to DataFrames, increasing mutation support (:pr:`1858`)

Delayed
^^^^^^^

- Changed behaviour for ``delayed(nout=0)`` and ``delayed(nout=1)``:
  ``delayed(nout=1)`` does not default to ``out=None`` anymore, and
  ``delayed(nout=0)`` is also enabled. I.e. functions with return
  tuples of length 1 or 0 can be handled correctly. This is especially
  handy, if functions with a variable amount of outputs are wrapped by
  ``delayed``. E.g. a trivial example:
  ``delayed(lambda *args: args, nout=len(vals))(*vals)``

Core
^^^^

- Refactor core byte ingest (:pr:`1768`), (:pr:`1774`)
- Improve import time (:pr:`1833`)


.. _v0.12.0 / 2016-11-03:

0.12.0 / 2016-11-03
-------------------

DataFrame
^^^^^^^^^
- Return a series when functions given to ``dataframe.map_partitions`` return
  scalars (:pr:`1515`)
- Fix type size inference for series (:pr:`1513`)
- ``dataframe.DataFrame.categorize`` no longer includes missing values
  in the ``categories``. This is for compatibility with a `pandas change <https://github.com/pydata/pandas/pull/10929>`_ (:pr:`1565`)
- Fix head parser error in ``dataframe.read_csv`` when some lines have quotes
  (:pr:`1495`)
- Add ``dataframe.reduction`` and ``series.reduction`` methods to apply generic
  row-wise reduction to dataframes and series (:pr:`1483`)
- Add ``dataframe.select_dtypes``, which mirrors the `pandas method <https://pandas.pydata.org/pandas-docs/version/0.18.1/generated/pandas.DataFrame.select_dtypes.html>`_ (:pr:`1556`)
- ``dataframe.read_hdf`` now supports reading ``Series`` (:pr:`1564`)
- Support Pandas 0.19.0 (:pr:`1540`)
- Implement ``select_dtypes`` (:pr:`1556`)
- String accessor works with indexes (:pr:`1561`)
- Add pipe method to dask.dataframe (:pr:`1567`)
- Add ``indicator`` keyword to merge (:pr:`1575`)
- Support Series in ``read_hdf`` (:pr:`1575`)
- Support Categories with missing values (:pr:`1578`)
- Support inplace operators like ``df.x += 1`` (:pr:`1585`)
- Str accessor passes through args and kwargs (:pr:`1621`)
- Improved groupby support for single-machine multiprocessing scheduler
  (:pr:`1625`)
- Tree reductions (:pr:`1663`)
- Pivot tables (:pr:`1665`)
- Add clip (:pr:`1667`), align (:pr:`1668`), combine_first (:pr:`1725`), and
  any/all (:pr:`1724`)
- Improved handling of divisions on dask-pandas merges (:pr:`1666`)
- Add ``groupby.aggregate`` method (:pr:`1678`)
- Add ``dd.read_table`` function (:pr:`1682`)
- Improve support for multi-level columns (:pr:`1697`) (:pr:`1712`)
- Support 2d indexing in ``loc`` (:pr:`1726`)
- Extend ``resample`` to include DataFrames (:pr:`1741`)
- Support dask.array ufuncs on dask.dataframe objects (:pr:`1669`)


Array
^^^^^
- Add information about how ``dask.array`` ``chunks`` argument work (:pr:`1504`)
- Fix field access with non-scalar fields in ``dask.array`` (:pr:`1484`)
- Add concatenate= keyword to atop to concatenate chunks of contracted dimensions
- Optimized slicing performance (:pr:`1539`) (:pr:`1731`)
- Extend ``atop`` with a ``concatenate=`` (:pr:`1609`) ``new_axes=``
  (:pr:`1612`) and ``adjust_chunks=`` (:pr:`1716`) keywords
- Add clip (:pr:`1610`) swapaxes (:pr:`1611`) round (:pr:`1708`) repeat
- Automatically align chunks in ``atop``-backed operations (:pr:`1644`)
- Cull dask.arrays on slicing (:pr:`1709`)

Bag
^^^
- Fix issue with callables in ``bag.from_sequence`` being interpreted as
  tasks (:pr:`1491`)
- Avoid non-lazy memory use in reductions (:pr:`1747`)

Administration
^^^^^^^^^^^^^^

- Added changelog (:pr:`1526`)
- Create new threadpool when operating from thread (:pr:`1487`)
- Unify example documentation pages into one (:pr:`1520`)
- Add versioneer for git-commit based versions (:pr:`1569`)
- Pass through node_attr and edge_attr keywords in dot visualization
  (:pr:`1614`)
- Add continuous testing for Windows with Appveyor (:pr:`1648`)
- Remove use of multiprocessing.Manager (:pr:`1653`)
- Add global optimizations keyword to compute (:pr:`1675`)
- Micro-optimize get_dependencies (:pr:`1722`)


.. _v0.11.0 / 2016-08-24:

0.11.0 / 2016-08-24
-------------------

Major Points
^^^^^^^^^^^^

DataFrames now enforce knowing full metadata (columns, dtypes) everywhere.
Previously we would operate in an ambiguous state when functions lost dtype
information (such as ``apply``).  Now all dataframes always know their dtypes
and raise errors asking for information if they are unable to infer (which
they usually can).  Some internal attributes like ``_pd`` and
``_pd_nonempty`` have been moved.

The internals of the distributed scheduler have been refactored to
transition tasks between explicit states.  This improves resilience,
reasoning about scheduling, plugin operation, and logging.  It also makes
the scheduler code easier to understand for newcomers.

Breaking Changes
^^^^^^^^^^^^^^^^

- The ``distributed.s3`` and ``distributed.hdfs`` namespaces are gone.  Use
  protocols in normal methods like ``read_text('s3://...'`` instead.
- ``Dask.array.reshape`` now errs in some cases where previously it would have
  create a very large number of tasks


.. _v0.10.2 / 2016-07-27:

0.10.2 / 2016-07-27
-------------------

- More Dataframe shuffles now work in distributed settings, ranging from
  setting-index to hash joins, to sorted joins and groupbys.
- Dask passes the full test suite when run when under in Python's
  optimized-OO mode.
- On-disk shuffles were found to produce wrong results in some
  highly-concurrent situations, especially on Windows.  This has been resolved
  by a fix to the partd library.
- Fixed a growth of open file descriptors that occurred under large data
  communications
- Support ports in the ``--bokeh-whitelist`` option ot dask-scheduler to better
  routing of web interface messages behind non-trivial network settings
- Some improvements to resilience to worker failure (though other known
  failures persist)
- You can now start an IPython kernel on any worker for improved debugging and
  analysis
- Improvements to ``dask.dataframe.read_hdf``, especially when reading from
  multiple files and docs


.. _v0.10.0 / 2016-06-13:

0.10.0 / 2016-06-13
-------------------

Major Changes
^^^^^^^^^^^^^

- This version drops support for Python 2.6
- Conda packages are built and served from conda-forge
- The ``dask.distributed`` executables have been renamed from dfoo to dask-foo.
  For example dscheduler is renamed to dask-scheduler
- Both Bag and DataFrame include a preliminary distributed shuffle.

Bag
^^^

- Add task-based shuffle for distributed groupbys
- Add accumulate for cumulative reductions

DataFrame
^^^^^^^^^

- Add a task-based shuffle suitable for distributed joins, groupby-applys, and
  set_index operations.  The single-machine shuffle remains untouched (and
  much more efficient.)
- Add support for new Pandas rolling API with improved communication
  performance on distributed systems.
- Add ``groupby.std/var``
- Pass through S3/HDFS storage options in ``read_csv``
- Improve categorical partitioning
- Add eval, info, isnull, notnull for dataframes

Distributed
^^^^^^^^^^^

- Rename executables like dscheduler to dask-scheduler
- Improve scheduler performance in the many-fast-tasks case (important for
  shuffling)
- Improve work stealing to be aware of expected function run-times and data
  sizes.  The drastically increases the breadth of algorithms that can be
  efficiently run on the distributed scheduler without significant user
  expertise.
- Support maximum buffer sizes in streaming queues
- Improve Windows support when using the Bokeh diagnostic web interface
- Support compression of very-large-bytestrings in protocol
- Support clean cancellation of submitted futures in Joblib interface

Other
^^^^^

- All dask-related projects (dask, distributed, s3fs, hdfs, partd) are now
  building conda packages on conda-forge.
- Change credential handling in s3fs to only pass around delegated credentials
  if explicitly given secret/key.  The default now is to rely on managed
  environments.  This can be changed back by explicitly providing a keyword
  argument.  Anonymous mode must be explicitly declared if desired.


.. _v0.9.0 / 2016-05-11:

0.9.0 / 2016-05-11
------------------

API Changes
^^^^^^^^^^^

- ``dask.do`` and ``dask.value`` have been renamed to ``dask.delayed``
- ``dask.bag.from_filenames`` has been renamed to ``dask.bag.read_text``
- All S3/HDFS data ingest functions like ``db.from_s3`` or
  ``distributed.s3.read_csv`` have been moved into the plain ``read_text``,
  ``read_csv functions``, which now support protocols, like
  ``dd.read_csv('s3://bucket/keys*.csv')``

Array
^^^^^

- Add support for ``scipy.LinearOperator``
- Improve optional locking to on-disk data structures
- Change rechunk to expose the intermediate chunks

Bag
^^^

- Rename ``from_filename``\ s to ``read_text``
- Remove ``from_s3`` in favor of ``read_text('s3://...')``

DataFrame
^^^^^^^^^

- Fixed numerical stability issue for correlation and covariance
- Allow no-hash ``from_pandas`` for speedy round-trips to and from-pandas
  objects
- Generally reengineered ``read_csv`` to be more in line with Pandas behavior
- Support fast ``set_index`` operations for sorted columns

Delayed
^^^^^^^

- Rename ``do/value`` to ``delayed``
- Rename ``to/from_imperative`` to ``to/from_delayed``

Distributed
^^^^^^^^^^^

- Move s3 and hdfs functionality into the dask repository
- Adaptively oversubscribe workers for very fast tasks
- Improve PyPy support
- Improve work stealing for unbalanced workers
- Scatter data efficiently with tree-scatters

Other
^^^^^

- Add lzma/xz compression support
- Raise a warning when trying to split unsplittable compression types, like
  gzip or bz2
- Improve hashing for single-machine shuffle operations
- Add new callback method for start state
- General performance tuning


.. _v0.8.1 / 2016-03-11:

0.8.1 / 2016-03-11
------------------

Array
^^^^^

- Bugfix for range slicing that could periodically lead to incorrect results.
- Improved support and resiliency of ``arg`` reductions (``argmin``, ``argmax``, etc.)

Bag
^^^

- Add ``zip`` function

DataFrame
^^^^^^^^^

- Add ``corr`` and ``cov`` functions
- Add ``melt`` function
- Bugfixes for io to bcolz and hdf5


.. _v0.8.0 / 2016-02-20:

0.8.0 / 2016-02-20
------------------

Array
^^^^^

- Changed default array reduction split from 32 to 4
- Linear algebra, ``tril``, ``triu``, ``LU``, ``inv``, ``cholesky``,
  ``solve``, ``solve_triangular``, ``eye``, ``lstsq``, ``diag``, ``corrcoef``.

Bag
^^^

- Add tree reductions
- Add range function
- drop ``from_hdfs`` function (better functionality now exists in hdfs3 and
  distributed projects)

DataFrame
^^^^^^^^^

- Refactor ``dask.dataframe`` to include a full empty pandas dataframe as
  metadata.  Drop the ``.columns`` attribute on Series
- Add Series categorical accessor, series.nunique, drop the ``.columns``
  attribute for series.
- ``read_csv`` fixes (multi-column parse_dates, integer column names, etc. )
- Internal changes to improve graph serialization

Other
^^^^^

- Documentation updates
- Add from_imperative and to_imperative functions for all collections
- Aesthetic changes to profiler plots
- Moved the dask project to a new dask organization


.. _v0.7.6 / 2016-01-05:

0.7.6 / 2016-01-05
------------------

Array
^^^^^
- Improve thread safety
- Tree reductions
- Add ``view``, ``compress``, ``hstack``, ``dstack``, ``vstack`` methods
- ``map_blocks`` can now remove and add dimensions

DataFrame
^^^^^^^^^
- Improve thread safety
- Extend sampling to include replacement options

Imperative
^^^^^^^^^^
- Removed optimization passes that fused results.

Core
^^^^

- Removed ``dask.distributed``
- Improved performance of blocked file reading
- Serialization improvements
- Test Python 3.5


.. _v0.7.4 / 2015-10-23:

0.7.4 / 2015-10-23
------------------

This was mostly a bugfix release. Some notable changes:

- Fix minor bugs associated with the release of numpy 1.10 and pandas 0.17
- Fixed a bug with random number generation that would cause repeated blocks
  due to the birthday paradox
- Use locks in ``dask.dataframe.read_hdf`` by default to avoid concurrency
  issues
- Change ``dask.get`` to point to ``dask.async.get_sync`` by default
- Allow visualization functions to accept general graphviz graph options like
  rankdir='LR'
- Add reshape and ravel to ``dask.array``
- Support the creation of ``dask.arrays`` from ``dask.imperative`` objects

Deprecation
^^^^^^^^^^^

This release also includes a deprecation warning for ``dask.distributed``, which
will be removed in the next version.

Future development in distributed computing for dask is happening here:
https://distributed.dask.org . General feedback on that project is most
welcome from this community.


.. _v0.7.3 / 2015-09-25:

0.7.3 / 2015-09-25
------------------

Diagnostics
^^^^^^^^^^^
- A utility for profiling memory and cpu usage has been added to the
  ``dask.diagnostics`` module.

DataFrame
^^^^^^^^^
This release improves coverage of the pandas API. Among other things
it includes ``nunique``, ``nlargest``, ``quantile``. Fixes encoding issues
with reading non-ascii csv files. Performance improvements and  bug fixes
with resample. More flexible read_hdf with globbing. And many more. Various
bug fixes in ``dask.imperative`` and ``dask.bag``.


.. _v0.7.0 / 2015-08-15:

0.7.0 / 2015-08-15
------------------

DataFrame
^^^^^^^^^
This release includes significant bugfixes and alignment with the Pandas API.
This has resulted both from use and from recent involvement by Pandas core
developers.

- New operations: query, rolling operations, drop
- Improved operations: quantiles, arithmetic on full dataframes, dropna,
  constructor logic, merge/join, elemwise operations, groupby aggregations

Bag
^^^
- Fixed a bug in fold where with a null default argument

Array
^^^^^
- New operations: da.fft module, da.image.imread

Infrastructure
^^^^^^^^^^^^^^
- The array and dataframe collections create graphs with deterministic keys.
  These tend to be longer (hash strings) but should be consistent between
  computations.  This will be useful for caching in the future.
- All collections (Array, Bag, DataFrame) inherit from common subclass


.. _v0.6.1 / 2015-07-23:

0.6.1 / 2015-07-23
------------------

Distributed
^^^^^^^^^^^
- Improved (though not yet sufficient) resiliency for ``dask.distributed``
  when workers die

DataFrame
^^^^^^^^^
- Improved writing to various formats, including to_hdf, to_castra, and
  to_csv
- Improved creation of dask DataFrames from dask Arrays and Bags
- Improved support for categoricals and various other methods

Array
^^^^^
- Various bug fixes
- Histogram function

Scheduling
^^^^^^^^^^
- Added tie-breaking ordering of tasks within parallel workloads to
  better handle and clear intermediate results

Other
^^^^^
- Added the dask.do function for explicit construction of graphs with
  normal python code
- Traded pydot for graphviz library for graph printing to support Python3
- There is also a gitter chat room and a stackoverflow tag


.. _`Guido Imperiale`: https://github.com/crusaderky
.. _`John A Kirkham`: https://github.com/jakirkham
.. _`Matthew Rocklin`: https://github.com/mrocklin
.. _`Jim Crist`: https://github.com/jcrist
.. _`James Bourbeau`: https://github.com/jrbourbeau
.. _`James Munroe`: https://github.com/jmunroe
.. _`Thomas Caswell`: https://github.com/tacaswell
.. _`Tom Augspurger`: https://github.com/tomaugspurger
.. _`Uwe Korn`: https://github.com/xhochy
.. _`Christopher Prohm`: https://github.com/chmp
.. _`@xwang777`: https://github.com/xwang777
.. _`@fjetter`: https://github.com/fjetter
.. _`@Ced4`: https://github.com/Ced4
.. _`Ian Hopkinson`: https://github.com/IanHopkinson
.. _`Stephan Hoyer`: https://github.com/shoyer
.. _`Albert DeFusco`: https://github.com/AlbertDeFusco
.. _`Markus Gonser`: https://github.com/magonser
.. _`Martijn Arts`: https://github.com/mfaafm
.. _`Jon Mease`: https://github.com/jonmmease
.. _`Xander Johnson`: https://github.com/metasyn
.. _`Nir`: https://github.com/nirizr
.. _`Keisuke Fujii`: https://github.com/fujiisoup
.. _`Roman Yurchak`: https://github.com/rth
.. _`Max Epstein`: https://github.com/MaxPowerWasTaken
.. _`Simon Perkins`: https://github.com/sjperkins
.. _`Richard Postelnik`: https://github.com/postelrich
.. _`Daniel Collins`: https://github.com/dancollins34
.. _`Gabriele Lanaro`: https://github.com/gabrielelanaro
.. _`Jörg Dietrich`: https://github.com/joergdietrich
.. _`Christopher Ren`: https://github.com/cr458
.. _`Martin Durant`: https://github.com/martindurant
.. _`Thrasibule`: https://github.com/thrasibule
.. _`Dieter Weber`: https://github.com/uellue
.. _`Apostolos Vlachopoulos`: https://github.com/avlahop
.. _`Jesse Vogt`: https://github.com/jessevogt
.. _`Pierre Bartet`: https://github.com/Pierre-Bartet
.. _`Scott Sievert`: https://github.com/stsievert
.. _`Jeremy Chen`: https://github.com/convexset
.. _`Marc Pfister`: https://github.com/drwelby
.. _`Matt Lee`: https://github.com/mathewlee11
.. _`Yu Feng`: https://github.com/rainwoodman
.. _`@andrethrill`: https://github.com/andrethrill
.. _`@beomi`: https://github.com/beomi
.. _`Henrique Ribeiro`: https://github.com/henriqueribeiro
.. _`Marco Rossi`: https://github.com/m-rossi
.. _`Itamar Turner-Trauring`: https://github.com/itamarst
.. _`Mike Neish`: https://github.com/neishm
.. _`Mark Harfouche`: https://github.com/hmaarrfk
.. _`George Sakkis`: https://github.com/gsakkis
.. _`Ziyao Wei`: https://github.com/ZiyaoWei
.. _`Jacob Tomlinson`: https://github.com/jacobtomlinson
.. _`Elliott Sales de Andrade`: https://github.com/QuLogic
.. _`Gerome Pistre`: https://github.com/GPistre
.. _`Cloves Almeida`: https://github.com/cjalmeida
.. _`Tobias de Jong`: https://github.com/tadejong
.. _`Irina Truong`: https://github.com/j-bennet
.. _`Eric Bonfadini`: https://github.com/eric-bonfadini
.. _`Danilo Horta`: https://github.com/horta
.. _`@hugovk`: https://github.com/hugovk
.. _`Jan Margeta`: https://github.com/jmargeta
.. _`John Mrziglod`: https://github.com/JohnMrziglod
.. _`Christoph Moehl`: https://github.com/cmohl2013
.. _`Anderson Banihirwe`: https://github.com/andersy005
.. _`Javad`: https://github.com/javad94
.. _`Daniel Rothenberg`: https://github.com/darothen
.. _`Hans Moritz Günther`: https://github.com/hamogu
.. _`@rtobar`: https://github.com/rtobar
.. _`Julia Signell`: https://github.com/jsignell
.. _`Sriharsha Hatwar`: https://github.com/Sriharsha-hatwar
.. _`Bruce Merry`: https://github.com/bmerry
.. _`Joe Hamman`: https://github.com/jhamman
.. _`Robert Sare`: https://github.com/rmsare
.. _`Jeremy Chan`: https://github.com/convexset
.. _`Eric Wolak`: https://github.com/epall
.. _`Miguel Farrajota`: https://github.com/farrajota
.. _`Zhenqing Li`: https://github.com/DigitalPig
.. _`Matthias Bussonier`: https://github.com/Carreau
.. _`Jan Koch`: https://github.com/datajanko
.. _`Bart Broere`: https://github.com/bartbroere
.. _`Rahul Vaidya`: https://github.com/rvaidya
.. _`Justin Dennison`: https://github.com/justin1dennison
.. _`Antonino Ingargiola`: https://github.com/tritemio
.. _`TakaakiFuruse`: https://github.com/TakaakiFuruse
.. _`samc0de`: https://github.com/samc0de
.. _`Armin Berres`: https://github.com/aberres
.. _`Damien Garaud`: https://github.com/geraud
.. _`Jonathan Fraine`: https://github.com/exowanderer
.. _`Carlos Valiente`: https://github.com/carletes
.. _`@milesial`: https://github.com/milesial
.. _`Paul Vecchio`: https://github.com/vecchp
.. _`Johnnie Gray`: https://github.com/jcmgray
.. _`Diane Trout`: https://github.com/detrout
.. _`Marco Neumann`: https://github.com/crepererum
.. _`Mina Farid`: https://github.com/minafarid
.. _`@slnguyen`: https://github.com/slnguyen
.. _`Gábor Lipták`: https://github.com/gliptak
.. _`David Hoese`: https://github.com/djhoese
.. _`Daniel Li`: https://github.com/li-dan
.. _`Prabakaran Kumaresshan`: https://github.com/nixphix
.. _`Daniel Saxton`: https://github.com/dsaxton
.. _`Jendrik Jördening`: https://github.com/jendrikjoe
.. _`Takahiro Kojima`: https://github.com/515hikaru
.. _`Stuart Berg`: https://github.com/stuarteberg
.. _`Guillaume Eynard-Bontemps`: https://github.com/guillaumeeb
.. _`Adam Beberg`: https://github.com/beberg
.. _`Roma Sokolov`: https://github.com/little-arhat
.. _`Daniel Severo`: https://github.com/dsevero
.. _`Michał Jastrzębski`: https://github.com/inc0
.. _`Janne Vuorela`: https://github.com/Dimplexion
.. _`Ross Petchler`: https://github.com/rpetchler
.. _`Aploium`: https://github.com/aploium
.. _`Peter Andreas Entschev`: https://github.com/pentschev
.. _`@JulianWgs`: https://github.com/JulianWgs
.. _`Shyam Saladi`: https://github.com/smsaladi
.. _`Joe Corbett`: https://github.com/jcorb
.. _`@HSR05`: https://github.com/HSR05
.. _`Benjamin Zaitlen`: https://github.com/quasiben
.. _`Brett Naul`: https://github.com/bnaul
.. _`Justin Poehnelt`: https://github.com/jpoehnelt
.. _`Dan O'Donovan`: https://github.com/danodonovan
.. _`amerkel2`: https://github.com/amerkel2
.. _`Justin Waugh`: https://github.com/bluecoconut
.. _`Brian Chu`: https://github.com/bchu
.. _`Álvaro Abella Bascarán`: https://github.com/alvaroabascar
.. _`Aaron Fowles`: https://github.com/aaronfowles
.. _`Søren Fuglede Jørgensen`: https://github.com/fuglede
.. _`Hameer Abbasi`: https://github.com/hameerabbasi
.. _`Philipp Rudiger`: https://github.com/philippjfr
.. _`gregrf`: https://github.com/gregrf
.. _`Ian Rose`: https://github.com/ian-r-rose
.. _`Genevieve Buckley`: https://github.com/GenevieveBuckley
.. _`Michael Eaton`: https://github.com/mpeaton
.. _`Isaiah Norton`: https://github.com/hnorton
.. _`Nick Becker`: https://github.com/beckernick
.. _`Nathan Matare`: https://github.com/nmatare
.. _`@asmith26`: https://github.com/asmith26
.. _`Abhinav Ralhan`: https://github.com/abhinavralhan
.. _`Christian Hudon`: https://github.com/chrish42
.. _`Alistair Miles`: https://github.com/alimanfoo
.. _`Henry Pinkard`: https://github.com/
.. _`Ian Bolliger`: https://github.com/bolliger32
.. _`Mark Bell`: https://github.com/MarkCBell
.. _`Cody Johnson`: https://github.com/codercody
.. _`Endre Mark Borza`: https://github.com/endremborza
.. _`asmith26`: https://github.com/asmith26
.. _`Philipp S. Sommer`: https://github.com/Chilipp
.. _`mcsoini`: https://github.com/mcsoini
.. _`Ksenia Bobrova`: https://github.com/almaleksia
.. _`tpanza`: https://github.com/tpanza
.. _`Richard J Zamora`: https://github.com/rjzamora
.. _`Lijo Jose`: https://github.com/lijose
.. _`btw08`: https://github.com/btw08
.. _`Jorge Pessoa`: https://github.com/jorge-pessoa
.. _`Guillaume Lemaitre`: https://github.com/glemaitre
.. _`Bouwe Andela`: https://github.com/bouweandela
.. _`mbarkhau`: https://github.com/mbarkhau
.. _`Hugo`: https://github.com/hugovk
.. _`Paweł Kordek`: https://github.com/kordek
.. _`Ralf Gommers`: https://github.com/rgommers
.. _`Davis Bennett`: https://github.com/d-v-b
.. _`Willi Rath`: https://github.com/willirath
.. _`David Brochart`: https://github.com/davidbrochart
.. _`GALI PREM SAGAR`: https://github.com/galipremsagar
.. _`tshatrov`: https://github.com/tshatrov
.. _`Dustin Tindall`: https://github.com/dustindall
.. _`Sean McKenna`: https://github.com/seanmck
.. _`msbrown47`: https://github.com/msbrown47
.. _`Natalya Rapstine`: https://github.com/natalya-patrikeeva
.. _`Loïc Estève`: https://github.com/lesteve
.. _`Xavier Holt`: https://github.com/xavi-ai
.. _`Sarah Bird`: https://github.com/birdsarah
.. _`Doug Davis`: https://github.com/douglasdavis
.. _`Nicolas Hug`: https://github.com/NicolasHug
.. _`Blane`: https://github.com/BlaneG
.. _`Ivars Geidans`: https://github.com/ivarsfg
.. _`Scott Sievert`: https://github.com/stsievert
.. _`estebanag`: https://github.com/estebanag
.. _`Benoit Bovy`: https://github.com/benbovy
.. _`Gabe Joseph`: https://github.com/gjoseph92
.. _`therhaag`: https://github.com/therhaag
.. _`Arpit Solanki`: https://github.com/arpit1997
.. _`Oliver Hofkens`: https://github.com/OliverHofkens
.. _`Hongjiu Zhang`: https://github.com/hongzmsft
.. _`Wes Roach`: https://github.com/WesRoach
.. _`DomHudson`: https://github.com/DomHudson
.. _`Eugene Huang`: https://github.com/eugeneh101
.. _`Christopher J. Wright`: https://github.com/CJ-Wright
.. _`Mahmut Bulut`: https://github.com/vertexclique
.. _`Ben Jeffery`: https://github.com/benjeffery
.. _`Ryan Nazareth`: https://github.com/ryankarlos
.. _`garanews`: https://github.com/garanews
.. _`Vijayant`: https://github.com/VijayantSoni
.. _`Ryan Abernathey`: https://github.com/rabernat
.. _`Norman Barker`: https://github.com/normanb
.. _`darindf`: https://github.com/darindf
.. _`Ryan Grout`: https://github.com/groutr
.. _`Krishan Bhasin`: https://github.com/KrishanBhasin
.. _`Albert DeFusco`: https://github.com/AlbertDeFusco
.. _`Bruno Bonfils`: https://github.com/asyd
.. _`Petio Petrov`: https://github.com/petioptrv
.. _`Mads R. B. Kristensen`: https://github.com/madsbk
.. _`Prithvi MK`: https://github.com/pmk21
.. _`Eric Dill`: https://github.com/ericdill
.. _`Gina Helfrich`: https://github.com/Dr-G
.. _`ossdev07`: https://github.com/ossdev07
.. _`Nuno Gomes Silva`: https://github.com/mgsnuno
.. _`Ray Bell`: https://github.com/raybellwaves
.. _`Deepak Cherian`: https://github.com/dcherian
.. _`Matteo De Wint`: https://github.com/mdwint
.. _`Tim Gates`: https://github.com/timgates42
.. _`Erik Welch`: https://github.com/eriknw
.. _`Christian Wesp`: https://github.com/ChrWesp
.. _`Shiva Raisinghani`: https://github.com/exemplary-citizen
.. _`Thomas A Caswell`: https://github.com/tacaswell
.. _`Timost`: https://github.com/Timost
.. _`Maarten Breddels`: https://github.com/maartenbreddels
.. _`Devin Petersohn`: https://github.com/devin-petersohn
.. _`dfonnegra`: https://github.com/dfonnegra
.. _`Chris Roat`: https://github.com/ChrisRoat
.. _`H. Thomson Comer`: https://github.com/thomcom
.. _`Gerrit Holl`: https://github.com/gerritholl
.. _`Thomas Robitaille`: https://github.com/astrofrog
.. _`Yifan Gu`: https://github.com/gyf304
.. _`Surya Avala`: https://github.com/suryaavala
.. _`Cyril Shcherbin`: https://github.com/shcherbin
.. _`Ram Rachum`: https://github.com/cool-RR
.. _`Igor Gotlibovych`: https://github.com/ig248
.. _`K.-Michael Aye`: https://github.com/michaelaye
.. _`Yetunde Dada`: https://github.com/yetudada
.. _`Andrew Thomas`: https://github.com/amcnicho
.. _`rockwellw`: https://github.com/rockwellw
.. _`Gil Forsyth`: https://github.com/gforsyth
.. _`Thomas J. Fan`: https://github.com/thomasjpfan
.. _`Henrik Andersson`: https://github.com/hnra
.. _`James Lamb`: https://github.com/jameslamb
.. _`Corey J. Nolet`: https://github.com/cjnolet
.. _`Chuanzhu Xu`: https://github.com/xcz011
.. _`Lucas Rademaker`: https://github.com/lr4d
.. _`JulianWgs`: https://github.com/JulianWgs
.. _`psimaj`: https://github.com/psimaj
.. _`mlondschien`: https://github.com/mlondschien
.. _`petiop`: https://github.com/petiop
.. _`Richard (Rick) Zamora`: https://github.com/rjzamora
.. _`Mark Boer`: https://github.com/mark-boer
.. _`Florian Jetter`: https://github.com/fjetter
.. _`Adam Lewis`: https://github.com/Adam-D-Lewis
.. _`David Chudzicki`: https://github.com/dchudz
.. _`Nick Evans`: https://github.com/nre
.. _`Kai Mühlbauer`: https://github.com/kmuehlbauer
.. _`swapna`: https://github.com/swapna-pg
.. _`Antonio Ercole De Luca`: https://github.com/eracle
.. _`Amol Umbarkar`: https://github.com/mindhash
.. _`noreentry`: https://github.com/noreentry
.. _`Marius van Niekerk`: https://github.com/mariusvniekerk
.. _`Tung Dang`: https://github.com/3cham
.. _`Jim Crist-Harif`: https://github.com/jcrist
.. _`Brian Larsen`: https://github.com/brl0
.. _`Nils Braun`: https://github.com/nils-braun
.. _`Scott Sanderson`: https://github.com/ssanderson
.. _`Gaurav Sheni`: https://github.com/gsheni
.. _`Andrew Fulton`: https://github.com/andrewfulton9
.. _`Stephanie Gott`: https://github.com/stephaniegott
.. _`Huite`: https://github.com/Huite
.. _`Ryan Williams`: https://github.com/ryan-williams
.. _`Eric Czech`: https://github.com/eric-czech
.. _`Abdulelah Bin Mahfoodh`: https://github.com/abduhbm
.. _`Ben Shaver`: https://github.com/bpshaver
.. _`Matthias Bussonnier`: https://github.com/Carreau
.. _`johnomotani`: https://github.com/johnomotani
.. _`Roberto Panai`: https://github.com/rpanai
.. _`Clark Zinzow`: https://github.com/clarkzinzow
.. _`Tom McTiernan`: https://github.com/tmct
.. _`joshreback`: https://github.com/joshreback
.. _`Jun Han (Johnson) Ooi`: https://github.com/tebesfinwo
.. _`Jim Circadian`: https://github.com/JimCircadian
.. _`Jack Xiaosong Xu`: https://github.com/jackxxu
.. _`Mike McCarty`: https://github.com/mmccarty
.. _`michaelnarodovitch`: https://github.com/michaelnarodovitch
.. _`David Sheldon`: https://github.com/davidsmf
.. _`McToel`: https://github.com/McToel
.. _`Kilian Lieret`: https://github.com/klieret
.. _`Noah D. Brenowitz`: https://github.com/nbren12
.. _`Jon Thielen`: https://github.com/jthielen
.. _`Poruri Sai Rahul`: https://github.com/rahulporuri
.. _`Kyle Nicholson`: https://github.com/kylejn27
.. _`Rafal Wojdyla`: https://github.com/ravwojdyla
.. _`Sam Grayson`: https://github.com/charmoniumQ
.. _`Madhur Tandon`: https://github.com/madhur-tandon
.. _`Joachim B Haga`: https://github.com/jobh
.. _`Pav A`: https://github.com/rs2
.. _`GFleishman`: https://github.com/GFleishman
.. _`Shang Wang`: https://github.com/shangw-nvidia
.. _`Illviljan`: https://github.com/Illviljan
.. _`Jan Borchmann`: https://github.com/jborchma
.. _`Ruben van de Geer`: https://github.com/rubenvdg
.. _`Akira Naruse`: https://github.com/anaruse
.. _`Zhengnan Zhao`: https://github.com/zzhengnan
.. _`Greg Hayes`: https://github.com/hayesgb
.. _`RogerMoens`: https://github.com/RogerMoens
.. _`manuels`: https://github.com/manuels
.. _`Rockwell Weiner`: https://github.com/rockwellw
.. _`Devanshu Desai`: https://github.com/devanshuDesai
.. _`David Katz`: https://github.com/DavidKatz-il
.. _`Stephannie Jimenez Gacha`: https://github.com/steff456
.. _`Magnus Nord`: https://github.com/magnunor
.. _`Callum Noble`: https://github.com/callumanoble
.. _`Pascal Bourgault`: https://github.com/aulemahal
.. _`Joris Van den Bossche`: https://github.com/jorisvandenbossche
.. _`Mark`: https://github.com/mchi
.. _`Kumar Bharath Prabhu`: https://github.com/kumarprabhu1988
.. _`Rob Malouf`: https://github.com/rmalouf
.. _`sdementen`: https://github.com/sdementen
.. _`patquem`: https://github.com/patquem
.. _`Amit Kumar`: https://github.com/aktech
.. _`D-Stacks`: https://github.com/D-Stacks
.. _`Kyle Barron`: https://github.com/kylebarron
.. _`Julius Busecke`: https://github.com/jbusecke
.. _`Sinclair Target`: https://github.com/sinclairtarget
.. _`Ashwin Srinath`: https://github.com/shwina
.. _`David Hassell`: https://github.com/davidhassell
.. _`brandon-b-miller`: https://github.com/brandon-b-miller
.. _`Hristo Georgiev`: https://github.com/hristog
.. _`Trevor Manz`: https://github.com/manzt
.. _`Madhu94`: https://github.com/Madhu94
.. _`gerrymanoim`: https://github.com/gerrymanoim
.. _`rs9w33`: https://github.com/rs9w33
.. _`Tom White`: https://github.com/tomwhite
.. _`Eoin Shanaghy`: https://github.com/eoinsha
.. _`Nick Vazquez`: https://github.com/nickvazz
.. _`cameron16`: https://github.com/cameron16
.. _`Daniel Mesejo-León`: https://github.com/mesejo
.. _`Naty Clementi`: https://github.com/ncclementi
.. _`JSKenyon`: https://github.com/jskenyon
.. _`Freyam Mehta`: https://github.com/freyam
.. _`Jiaming Yuan`: https://github.com/trivialfis
.. _`c-thiel`: https://github.com/c-thiel
.. _`Andrew Champion`: https://github.com/aschampion
.. _`keewis`: https://github.com/keewis
.. _`Maisie Marshall`: https://github.com/maisiemarshall
.. _`Vibhu Jawa`: https://github.com/VibhuJawa
.. _`Boaz Mohar`: https://github.com/boazmohar
.. _`Kristopher Overholt`: https://github.com/koverholt
.. _`tsuga`: https://github.com/tsuga
.. _`Gabriel Miretti`: https://github.com/gmiretti
.. _`Geoffrey Lentner`: https://github.com/glentner
.. _`Charles Blackmon-Luca`: https://github.com/charlesbluca
.. _`Bryan Van de Ven`: https://github.com/bryevdv
.. _`Fabian Gebhart`: https://github.com/fgebhart
.. _`Ross`: https://github.com/rhjmoore
.. _`gurunath`: https://github.com/rajagurunath
.. _`aa1371`: https://github.com/aa1371
.. _`Gregory R. Lee`: https://github.com/grlee77
.. _`Louis Maddox`: https://github.com/lmmx
.. _`Dahn`: https://github.com/DahnJ
.. _`Jordan Jensen`: https://github.com/dotNomad
.. _`Martin Fleischmann`: https://github.com/martinfleis
.. _`Robert Hales`: https://github.com/robalar
.. _`João Paulo Lacerda`: https://github.com/jopasdev
.. _`SnkSynthesis`: https://github.com/SnkSynthesis
.. _`JoranDox`: https://github.com/JoranDox
.. _`Kinshuk Dua`: https://github.com/kinshukdua
.. _`Suriya Senthilkumar`: https://github.com/suriya-it19
.. _`Vũ Trung Đức`: https://github.com/vutrungduc7593
.. _`Nathan Danielsen`: https://github.com/ndanielsen
.. _`Wallace Reis`: https://github.com/wreis
.. _`German Shiklov`: https://github.com/Jeremaiha-xmetix
.. _`Pankaj Patil`: https://github.com/Patil2099
.. _`Samuel Gaist`: https://github.com/sgaist
.. _`Marcel Coetzee`: https://github.com/marcelned
.. _`Matthew Powers`: https://github.com/MrPowers
.. _`Vyas Ramasubramani`: https://github.com/vyasr
.. _`Ayush Dattagupta`: https://github.com/ayushdg
.. _`FredericOdermatt`: https://github.com/FredericOdermatt
.. _`mihir`: https://github.com/ek234
.. _`Sarah Charlotte Johnson`: https://github.com/scharlottej13
.. _`ofirr`: https://github.com/ofirr
.. _`kori73`: https://github.com/kori73
.. _`TnTo`: https://github.com/TnTo
.. _`ParticularMiner`: https://github.com/ParticularMiner
.. _`aeisenbarth`: https://github.com/aeisenbarth
.. _`Aneesh Nema`: https://github.com/aneeshnema
.. _`Deepyaman Datta`: https://github.com/deepyaman
.. _`Maren Westermann`: https://github.com/marenwestermann
.. _`Michael Delgado`: https://github.com/delgadom
.. _`abergou`: https://github.com/abergou
.. _`Pavithra Eswaramoorthy`: https://github.com/pavithraes
.. _`Maxim Lippeveld`: https://github.com/MaximLippeveld
.. _`Kirito1397`: https://github.com/Kirito1397
API Reference
=============

Dask APIs generally follow from upstream APIs:

-  :doc:`Arrays<array-api>` follows NumPy
-  :doc:`DataFrames <dataframe-api>` follows Pandas
-  :doc:`Bag <bag-api>` follows map/filter/groupby/reduce common in Spark and Python iterators
-  :doc:`Delayed <delayed-api>` wraps general Python code
-  :doc:`Futures <futures>` follows `concurrent.futures <https://docs.python.org/3/library/concurrent.futures.html>`_ from the standard library for real-time computation.

.. toctree::
   :maxdepth: 1
   :hidden:

   Array <array-api.rst>
   DataFrame <dataframe-api.rst>
   Bag <bag-api.rst>
   Delayed <delayed-api.rst>
   Futures <futures>


Additionally, Dask has its own functions to start computations, persist data in
memory, check progress, and so forth that complement the APIs above.
These more general Dask functions are described below:

.. currentmodule:: dask

.. autosummary::
   compute
   is_dask_collection
   optimize
   persist
   visualize

These functions work with any scheduler.  More advanced operations are
available when using the newer scheduler and starting a
:obj:`dask.distributed.Client` (which, despite its name, runs nicely on a
single machine).  This API provides the ability to submit, cancel, and track
work asynchronously, and includes many functions for complex inter-task
workflows.  These are not necessary for normal operation, but can be useful for
real-time or advanced operation.

This more advanced API is available in the `Dask distributed documentation
<https://distributed.dask.org/en/latest/api.html>`_

.. autofunction:: annotate
.. autofunction:: compute
.. autofunction:: is_dask_collection
.. autofunction:: optimize
.. autofunction:: persist
.. autofunction:: visualize

Datasets
--------

Dask has a few helpers for generating demo datasets

.. currentmodule:: dask.datasets

.. autofunction:: make_people
.. autofunction:: timeseries

.. _api.utilities:

Utilities
---------

Dask has some public utility methods. These are primarily used for parsing
configuration values.

.. currentmodule:: dask.utils

.. autofunction:: format_bytes
.. autofunction:: format_time
.. autofunction:: parse_bytes
.. autofunction:: parse_timedelta
.. _array.chunks:

Chunks
======

Dask arrays are composed of many NumPy (or NumPy-like) arrays. How these arrays
are arranged can significantly affect performance.  For example, for a square
array you might arrange your chunks along rows, along columns, or in a more
square-like fashion. Different arrangements of NumPy arrays will be faster or
slower for different algorithms.

Thinking about and controlling chunking is important to optimize advanced
algorithms.

Specifying Chunk shapes
-----------------------

We always specify a ``chunks`` argument to tell dask.array how to break up the
underlying array into chunks.  We can specify ``chunks`` in a variety of ways:

1.  A uniform dimension size like ``1000``, meaning chunks of size ``1000`` in each dimension
2.  A uniform chunk shape like ``(1000, 2000, 3000)``, meaning chunks of size ``1000`` in the
    first axis, ``2000`` in the second axis, and ``3000`` in the third
3.  Fully explicit sizes of all blocks along all dimensions,
    like ``((1000, 1000, 500), (400, 400), (5, 5, 5, 5, 5))``
4.  A dictionary specifying chunk size per dimension like ``{0: 1000, 1: 2000,
    2: 3000}``.  This is just another way of writing the forms 2 and 3 above

Your chunks input will be normalized and stored in the third and most explicit
form.  Note that ``chunks`` stands for "chunk shape" rather than "number of
chunks", so specifying ``chunks=1`` means that you will have many chunks,
each with exactly one element.

For performance, a good choice of ``chunks`` follows the following rules:

1.  A chunk should be small enough to fit comfortably in memory.  We'll
    have many chunks in memory at once
2.  A chunk must be large enough so that computations on that chunk take
    significantly longer than the 1ms overhead per task that Dask scheduling
    incurs.  A task should take longer than 100ms
3.  Chunk sizes between 10MB-1GB are common, depending on the availability of
    RAM and the duration of computations
4.  Chunks should align with the computation that you want to do.

    For example, if you plan to frequently slice along a particular dimension,
    then it's more efficient if your chunks are aligned so that you have to
    touch fewer chunks.  If you want to add two arrays, then its convenient if
    those arrays have matching chunks patterns

5.  Chunks should align with your storage, if applicable.

    Array data formats are often chunked as well.  When loading or saving data,
    if is useful to have Dask array chunks that are aligned with the chunking
    of your storage, often an even multiple times larger in each direction


Unknown Chunks
--------------

Some arrays have unknown chunk sizes.  This arises whenever the size of an
array depends on lazy computations that we haven't yet performed like the
following:

.. code-block:: python

   >>> x = da.from_array(np.random.randn(100), chunks=20)
   >>> x += 0.1
   >>> y = x[x > 0]  # don't know how many values are greater than 0 ahead of time

Operations like the above result in arrays with unknown shapes and unknown
chunk sizes.  Unknown values within shape or chunks are designated using
``np.nan`` rather than an integer.  These arrays support many (but not all)
operations.  In particular, operations like slicing are not possible and will
result in an error.

.. code-block:: python

   >>> y.shape
   (np.nan,)
   >>> y[4]
   ...
   ValueError: Array chunk sizes unknown

   A possible solution: https://docs.dask.org/en/latest/array-chunks.html#unknown-chunks.
   Summary: to compute chunks sizes, use

       x.compute_chunk_sizes()  # for Dask Array
       ddf.to_dask_array(lengths=True)  # for Dask DataFrame ddf

Using :func:`~dask.array.Array.compute_chunk_sizes`  allows this example run:

.. code-block:: python

   >>> y.compute_chunk_sizes()
   dask.array<..., chunksize=(19,), ...>
   >>> y.shape
   (44,)
   >>> y[4].compute()
   0.78621774046566

Note that :func:`~dask.array.Array.compute_chunk_sizes` immediately performs computation and
modifies the array in-place.

Unknown chunksizes also occur when using a Dask DataFrame to create a Dask array:

.. code-block:: python

   >>> ddf = dask.dataframe.from_pandas(...)
   >>> ddf.to_dask_array()
   dask.array<..., shape=(nan, 2), ..., chunksize=(nan, 2)>

Using :func:`~dask.dataframe.DataFrame.to_dask_array` resolves this issue:

.. code-block:: python

   >>> ddf.to_dask_array(lengths=True)
   dask.array<..., shape=(100, 2), ..., chunksize=(20, 2)>

More details on :func:`~dask.dataframe.DataFrame.to_dask_array` are in mentioned in how to create a Dask
array from a Dask DataFrame in the :doc:`documentation on Dask array creation
<array-creation>`.

Chunks Examples
---------------

In this example we show how different inputs for ``chunks=`` cut up the following array::

   1 2 3 4 5 6
   7 8 9 0 1 2
   3 4 5 6 7 8
   9 0 1 2 3 4
   5 6 7 8 9 0
   1 2 3 4 5 6

Here, we show how different ``chunks=`` arguments split the array into different blocks

**chunks=3**: Symmetric blocks of size 3::

   1 2 3  4 5 6
   7 8 9  0 1 2
   3 4 5  6 7 8

   9 0 1  2 3 4
   5 6 7  8 9 0
   1 2 3  4 5 6

**chunks=2**: Symmetric blocks of size 2::

   1 2  3 4  5 6
   7 8  9 0  1 2

   3 4  5 6  7 8
   9 0  1 2  3 4

   5 6  7 8  9 0
   1 2  3 4  5 6

**chunks=(3, 2)**: Asymmetric but repeated blocks of size ``(3, 2)``::

   1 2  3 4  5 6
   7 8  9 0  1 2
   3 4  5 6  7 8

   9 0  1 2  3 4
   5 6  7 8  9 0
   1 2  3 4  5 6

**chunks=(1, 6)**: Asymmetric but repeated blocks of size ``(1, 6)``::

   1 2 3 4 5 6

   7 8 9 0 1 2

   3 4 5 6 7 8

   9 0 1 2 3 4

   5 6 7 8 9 0

   1 2 3 4 5 6

**chunks=((2, 4), (3, 3))**: Asymmetric and non-repeated blocks::

   1 2 3  4 5 6
   7 8 9  0 1 2

   3 4 5  6 7 8
   9 0 1  2 3 4
   5 6 7  8 9 0
   1 2 3  4 5 6

**chunks=((2, 2, 1, 1), (3, 2, 1))**: Asymmetric and non-repeated blocks::

   1 2 3  4 5  6
   7 8 9  0 1  2

   3 4 5  6 7  8
   9 0 1  2 3  4

   5 6 7  8 9  0

   1 2 3  4 5  6

**Discussion**

The latter examples are rarely provided by users on original data but arise from complex slicing and broadcasting operations.  Generally people use the simplest form until they need more complex forms.  The choice of chunks should align with the computations you want to do.

For example, if you plan to take out thin slices along the first dimension, then you might want to make that dimension skinnier than the others.  If you plan to do linear algebra, then you might want more symmetric blocks.


Loading Chunked Data
--------------------

Modern NDArray storage formats like HDF5, NetCDF, TIFF, and Zarr, allow arrays
to be stored in chunks or tiles so that blocks of data can be pulled out
efficiently without having to seek through a linear data stream.  It is best to
align the chunks of your Dask array with the chunks of your underlying data
store.

However, data stores often chunk more finely than is ideal for Dask array, so
it is common to choose a chunking that is a multiple of your storage chunk
size, otherwise you might incur high overhead.

For example, if you are loading a data store that is chunked in blocks of
``(100, 100)``, then you might choose a chunking more like ``(1000, 2000)`` that
is larger, but still evenly divisible by ``(100, 100)``.  Data storage
technologies will be able to tell you how their data is chunked.


Rechunking
----------

.. currentmodule:: dask.array

.. autosummary:: rechunk

Sometimes you need to change the chunking layout of your data.  For example,
perhaps it comes to you chunked row-wise, but you need to do an operation that
is much faster if done across columns.  You can change the chunking with the
``rechunk`` method.

.. code-block:: python

   x = x.rechunk((50, 1000))

Rechunking across axes can be expensive and incur a lot of communication, but
Dask array has fairly efficient algorithms to accomplish this.

You can pass rechunk any valid chunking form:

.. code-block:: python

   x = x.rechunk(1000)
   x = x.rechunk((50, 1000))
   x = x.rechunk({0: 50, 1: 1000})


.. _array-chunks.reshaping:

Reshaping
---------

The efficiency of :func:`dask.array.reshape` can depend strongly on the chunking
of the input array. In reshaping operations, there's the concept of "fast-moving"
or "high" axes. For a 2d array the second axis (``axis=1``) is the fastest-moving,
followed by the first. This means that if we draw a line indicating how values
are filled, we move across the "columns" first (along ``axis=1``), and then down
to the next row. Consider ``np.ones((3, 4)).reshape(12)``:

.. image:: images/reshape.png
   :alt: Visual representation of a 2-dimensional (3 rows by 4 colurmns) NumPy array being reshaped to 1 dimension (12 columns by 1 row). Arrows indicate the order in which values from the original array are copied to the new array, moving across the columns in axis 1 first before moving down to the next row in axis 0.

Now consider the impact of Dask's chunking on this operation. If the slow-moving
axis (just ``axis=0`` in this case) has chunks larger than size 1, we run into
a problem.

.. image:: images/reshape_problem.png

The first block has a shape ``(2, 2)``. Following the rules of ``reshape`` we
take the two values from the first row of block 1. But then we cross a chunk
boundary (from 1 to 2) while we still have two "unused" values in the first
block. There's no way to line up the input blocks with the output shape. We
need to somehow rechunk the input to be compatible with the output shape. We
have two options

1. Merge chunks using the logic in :meth:`dask.array.rechunk`. This avoids
   making two many tasks / blocks, at the cost of some communication and
   larger intermediates. This is the default behavior.
2. Use ``da.reshape(x, shape, merge_chunks=False)`` to avoid merging chunks
   by *splitting the input*. In particular, we can rechunk all the
   slow-moving axes to have a chunksize of 1. This avoids
   communication and moving around large amounts of data, at the cost of
   a larger task graph (potentially much larger, since the number of chunks
   on the slow-moving axes will equal the length of those axes.).

Visually, here's the second option:

.. image:: images/reshape_rechunked.png

Which if these is better depends on your problem. If communication is very
expensive and your data is relatively small along the slow-moving axes, then
``merge_chunks=False`` may be better. Let's compare the task graphs of these
two on a problem reshaping a 3-d array to a 2-d, where the input array doesn't
have ``chunksize=1`` on the slow-moving axes.

.. code-block:: python

   >>> a = da.from_array(np.arange(24).reshape(2, 3, 4), chunks=((2,), (2, 1), (2, 2)))
   >>> a
   dask.array<array, shape=(2, 3, 4), dtype=int64, chunksize=(2, 2, 2), chunktype=numpy.ndarray>
   >>> a.reshape(6, 4).visualize()

.. image:: images/merge_chunks.png

.. code-block:: python

   >>> a.reshape(6, 4, merge_chunks=False).visualize()

.. image:: images/merge_chunks_false.png

By default, some intermediate chunks chunks are merged, leading to a more complicated task
graph. With ``merge_chunks=False`` we split the input chunks (leading to more overall tasks,
depending on the size of the array) but avoid later communication.

Automatic Chunking
------------------

Chunks also includes three special values:

1.  ``-1``: no chunking along this dimension
2.  ``None``: no change to the chunking along this dimension (useful for rechunk)
3.  ``"auto"``: allow the chunking in this dimension to accommodate ideal chunk sizes

So, for example, one could rechunk a 3D array to have no chunking along the zeroth
dimension, but still have sensible chunk sizes as follows:

.. code-block:: python

   x = x.rechunk({0: -1, 1: 'auto', 2: 'auto'})

Or one can allow *all* dimensions to be auto-scaled to get to a good chunk
size:

.. code-block:: python

   x = x.rechunk('auto')

Automatic chunking expands or contracts all dimensions marked with ``"auto"``
to try to reach chunk sizes with a number of bytes equal to the config value
``array.chunk-size``, which is set to 128MiB by default, but which you can
change in your :doc:`configuration <configuration>`.

.. code-block:: python

   >>> dask.config.get('array.chunk-size')
   '128MiB'

Automatic rechunking tries to respect the median chunk shape of the
auto-rescaled dimensions, but will modify this to accommodate the shape of the
full array (can't have larger chunks than the array itself) and to find
chunk shapes that nicely divide the shape.

These values can also be used when creating arrays with operations like
``dask.array.ones`` or ``dask.array.from_array``

.. code-block:: python

   >>> dask.array.ones((10000, 10000), chunks=(-1, 'auto'))
   dask.array<wrapped, shape=(10000, 10000), dtype=float64, chunksize=(10000, 1250), chunktype=numpy.ndarray>
Customize initialization
========================

Often we want to run custom code when we start up or tear down a scheduler or
worker.  We might do this manually with functions like ``Client.run`` or
``Client.run_on_scheduler``, but this is error prone and difficult to automate.

To resolve this, Dask includes a few mechanisms to run arbitrary code around
the lifecycle of a Scheduler or Worker.

Preload Scripts
---------------

Both ``dask-scheduler`` and ``dask-worker`` support a ``--preload`` option that
allows custom initialization of each scheduler/worker respectively. A module or
Python file passed as a ``--preload`` value is guaranteed to be imported before
establishing any connection. A ``dask_setup(service)`` function is called if
found, with a ``Scheduler`` or ``Worker`` instance as the argument. As the
service stops, ``dask_teardown(service)`` is called if present.

To support additional configuration, a single ``--preload`` module may register
additional command-line arguments by exposing ``dask_setup`` as a  Click_
command.  This command will be used to parse additional arguments provided to
``dask-worker`` or ``dask-scheduler`` and will be called before service
initialization.

.. _Click: http://click.pocoo.org/

Example
~~~~~~~

As an example, consider the following file that creates a
`scheduler plugin <https://distributed.dask.org/en/latest/plugins.html>`_
and registers it with the scheduler

.. code-block:: python

   # scheduler-setup.py
   import click

   from distributed.diagnostics.plugin import SchedulerPlugin

   class MyPlugin(SchedulerPlugin):
       def __init__(self, print_count):
         self.print_count = print_count
         super().__init__()

       def add_worker(self, scheduler=None, worker=None, **kwargs):
           print("Added a new worker at:", worker)
           if self.print_count and scheduler is not None:
               print("Total workers:", len(scheduler.workers))

   @click.command()
   @click.option("--print-count/--no-print-count", default=False)
   def dask_setup(scheduler, print_count):
       plugin = MyPlugin(print_count)
       scheduler.add_plugin(plugin)

We can then run this preload script by referring to its filename (or module name
if it is on the path) when we start the scheduler::

   dask-scheduler --preload scheduler-setup.py --print-count

Types
~~~~~

Preloads can be specified as any of the following forms:

-   A path to a script, like ``/path/to/myfile.py``
-   A module name that is on the path, like ``my_module.initialize``
-   The text of a Python script, like ``import os; os.environ["A"] = "value"``

Configuration
~~~~~~~~~~~~~

Preloads can also be registered with configuration at the following values:

.. code-block:: yaml

   distributed:
     scheduler:
       preload:
       - "import os; os.environ['A'] = 'b'"  # use Python text
       - /path/to/myfile.py                  # or a filename
       - my_module                           # or a module name
       preload_argv:
       - []                                  # Pass optional keywords
       - ["--option", "value"]
       - []
     worker:
       preload: []
       preload_argv: []
     nanny:
       preload: []
       preload_argv: []

.. note::

   Because the ``dask-worker`` command needs to accept keywords for both the
   Worker and the Nanny (if a nanny is used) it has both a ``--preload`` and
   ``--preload-nanny`` keyword.  All extra keywords (like ``--print-count``
   above) will be sent to the workers rather than the nanny.  There is no way
   to specify extra keywords to the nanny preload scripts on the command line.
   We recommend the use of the more flexible configuration if this is
   necessary.


Worker Lifecycle Plugins
------------------------

You can also create a class with ``setup``, ``teardown``, and ``transition`` methods,
and register that class with the scheduler to give to every worker using the
``Client.register_worker_plugin`` method.

.. currentmodule:: distributed

.. autosummary::
   Client.register_worker_plugin

.. automethod:: Client.register_worker_plugin
   :noindex:
Debug
=====

Debugging parallel programs is hard.  Normal debugging tools like logging and
using ``pdb`` to interact with tracebacks stop working normally when exceptions
occur in far-away machines, different processes, or threads.

Dask has a variety of mechanisms to make this process easier.  Depending on
your situation, some of these approaches may be more appropriate than others.

These approaches are ordered from lightweight or easy solutions to more
involved solutions.

Exceptions
----------

When a task in your computation fails, the standard way of understanding what
went wrong is to look at the exception and traceback.  Often people do this
with the ``pdb`` module, IPython ``%debug`` or ``%pdb`` magics, or by just
looking at the traceback and investigating where in their code the exception
occurred.

Normally when a computation executes in a separate thread or a different
machine, these approaches break down.  To address this, Dask provides a few
mechanisms to recreate the normal Python debugging experience.

Inspect Exceptions and Tracebacks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, Dask already copies the exception and traceback wherever they
occur and reraises that exception locally.  If your task failed with a
``ZeroDivisionError`` remotely, then you'll get a ``ZeroDivisionError`` in your
interactive session.  Similarly you'll see a full traceback of where this error
occurred, which, just like in normal Python, can help you to identify the
troublesome spot in your code.

However, you cannot use the ``pdb`` module or ``%debug`` IPython magics with
these tracebacks to look at the value of variables during failure.  You can
only inspect things visually.  Additionally, the top of the traceback may be
filled with functions that are Dask-specific and not relevant to your
problem, so you can safely ignore these.

Both the single-machine and distributed schedulers do this.


Use the Single-Threaded Scheduler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dask ships with a simple single-threaded scheduler.  This doesn't offer any
parallel performance improvements but does run your Dask computation
faithfully in your local thread, allowing you to use normal tools like ``pdb``,
``%debug`` IPython magics, the profiling tools like the ``cProfile`` module, and
`snakeviz <https://jiffyclub.github.io/snakeviz/>`_.  This allows you to use
all of your normal Python debugging tricks in Dask computations, as long as you
don't need parallelism.

The single-threaded scheduler can be used, for example, by setting
``scheduler='single-threaded'`` in a compute call:

.. code-block:: python

    >>> x.compute(scheduler='single-threaded')

For more ways to configure schedulers, see the :ref:`scheduler configuration
documentation <scheduling-configuration>`.

This only works for single-machine schedulers.  It does not work with
``dask.distributed`` unless you are comfortable using the Tornado API (look at the
`testing infrastructure
<https://distributed.dask.org/en/latest/develop.html#writing-tests>`_
docs, which accomplish this).  Also, because this operates on a single machine,
it assumes that your computation can run on a single machine without exceeding
memory limits.  It may be wise to use this approach on smaller versions of your
problem if possible.


Rerun Failed Task Locally
~~~~~~~~~~~~~~~~~~~~~~~~~

If a remote task fails, we can collect the function and all inputs, bring them
to the local thread, and then rerun the function in hopes of triggering the
same exception locally where normal debugging tools can be used.

With the single-machine schedulers, use the ``rerun_exceptions_locally=True``
keyword:

.. code-block:: python

   >>> x.compute(rerun_exceptions_locally=True)

On the distributed scheduler use the ``recreate_error_locally`` method on
anything that contains ``Futures``:

.. code-block:: python

   >>> x.compute()
   ZeroDivisionError(...)

   >>> %pdb
   >>> future = client.compute(x)
   >>> client.recreate_error_locally(future)


Remove Failed Futures Manually
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes only parts of your computations fail, for example, if some rows of a
CSV dataset are faulty in some way.  When running with the distributed
scheduler, you can remove chunks of your data that have produced bad results if
you switch to dealing with Futures:

.. code-block:: python

   >>> import dask.dataframe as dd
   >>> df = ...           # create dataframe
   >>> df = df.persist()  # start computing on the cluster

   >>> from distributed.client import futures_of
   >>> futures = futures_of(df)  # get futures behind dataframe
   >>> futures
   [<Future: status: finished, type: pd.DataFrame, key: load-1>
    <Future: status: finished, type: pd.DataFrame, key: load-2>
    <Future: status: error, key: load-3>
    <Future: status: pending, key: load-4>
    <Future: status: error, key: load-5>]

   >>> # wait until computation is done
   >>> while any(f.status == 'pending' for f in futures):
   ...     sleep(0.1)

   >>> # pick out only the successful futures and reconstruct the dataframe
   >>> good_futures = [f for f in futures if f.status == 'finished']
   >>> df = dd.from_delayed(good_futures, meta=df._meta)

This is a bit of a hack, but often practical when first exploring messy data.
If you are using the concurrent.futures API (map, submit, gather), then this
approach is more natural.


Inspect Scheduling State
------------------------

Not all errors present themselves as exceptions.  For example, in a distributed
system workers may die unexpectedly, your computation may be unreasonably
slow due to inter-worker communication or scheduler overhead, or one of several
other issues.  Getting feedback about what's going on can help to identify
both failures and general performance bottlenecks.

For the single-machine scheduler, see :doc:`local diagnostic
<../diagnostics-local>` documentation.  The rest of the section will
assume that you are using the `distributed scheduler
<https://distributed.dask.org/en/latest/>`_ where these issues arise more
commonly.

Web Diagnostics
~~~~~~~~~~~~~~~

First, the distributed scheduler has a number of `diagnostic tools
<https://distributed.dask.org/en/latest/diagnosing-performance.html>`_ showing dozens of
recorded metrics like CPU, memory, network, and disk use, a history of previous
tasks, allocation of tasks to workers, worker memory pressure, work stealing,
open file handle limits, etc.  *Many* problems can be correctly diagnosed by
inspecting these pages.  By default, these are available at
``http://scheduler:8787/``, ``http://scheduler:8788/``, and ``http://worker:8789/``,
where ``scheduler`` and ``worker`` should be replaced by the addresses of the
scheduler and each of the workers. See `diagnosing performance docs
<https://distributed.dask.org/en/latest/diagnosing-performance.html>`_ for more information.

Logs
~~~~

The scheduler, workers, and client all emits logs using `Python's standard
logging module <https://docs.python.org/3/library/logging.html>`_.  By default,
these emit to standard error.  When Dask is launched by a cluster job scheduler
(SGE/SLURM/YARN/Mesos/Marathon/Kubernetes/whatever), that system will track
these logs and will have an interface to help you access them.  If you are
launching Dask on your own, they will probably dump to the screen unless you
`redirect stderr to a file
<https://en.wikipedia.org/wiki/Redirection_(computing)#Redirecting_to_and_from_the_standard_file_handles>`_
.

You can control the logging verbosity in the :doc:`../configuration`, for example,
the ``~/.config/dask/*.yaml`` files.
Defaults currently look like the following:

.. code-block:: yaml

   logging:
     distributed: info
     distributed.client: warning
     bokeh: error

Logging for specific components like ``distributed.client``,  ``distributed.scheduler``,
``distributed.nanny``,  ``distributed.worker``, etc. can each be independently configured.
So, for example, you could add a line like ``distributed.worker: debug`` to get
*very* verbose output from the workers.

Furthermore, you can explicitly assign handlers to loggers. The following example
assigns both file ("output.log") and console output to the scheduler and workers.
See the `python logging`_ documentation for information on the meaning of
specific terms here.

.. code-block:: yaml

    logging:
      version: 1
      handlers:
        file:
          class: logging.handlers.RotatingFileHandler
          filename: output.log
          level: INFO
        console:
          class: logging.StreamHandler
          level: INFO
      loggers:
        distributed.worker:
          level: INFO
          handlers:
            - file
            - console
        distributed.scheduler:
          level: INFO
          handlers:
            - file
            - console

.. _python logging: https://docs.python.org/3/library/logging.html


LocalCluster
------------

If you are using the distributed scheduler from a single machine, you may be
setting up workers manually using the command line interface or you may be
using `LocalCluster <https://distributed.dask.org/en/latest/api.html#cluster>`_
which is what runs when you just call ``Client()``:

.. code-block:: python

   >>> from dask.distributed import Client, LocalCluster
   >>> client = Client()  # This is actually the following two commands

   >>> cluster = LocalCluster()
   >>> client = Client(cluster.scheduler.address)

LocalCluster is useful because the scheduler and workers are in the same
process with you, so you can easily inspect their `state
<https://distributed.dask.org/en/latest/scheduling-state.html>`_ while
they run (they are running in a separate thread):

.. code-block:: python

   >>> cluster.scheduler.processing
   {'worker-one:59858': {'inc-123', 'add-443'},
    'worker-two:48248': {'inc-456'}}

You can also do this for the workers *if* you run them without nanny processes:

.. code-block:: python

   >>> cluster = LocalCluster(nanny=False)
   >>> client = Client(cluster)

This can be very helpful if you want to use the Dask distributed API and still
want to investigate what is going on directly within the workers.  Information
is not distilled for you like it is in the web diagnostics, but you have full
low-level access.


Inspect state with IPython
--------------------------

Sometimes you want to inspect the state of your cluster but you don't have the
luxury of operating on a single machine.  In these cases you can launch an
IPython kernel on the scheduler and on every worker, which lets you inspect
state on the scheduler and workers as computations are completing.

This does not give you the ability to run ``%pdb`` or ``%debug`` on remote
machines. The tasks are still running in separate threads, and so are not
easily accessible from an interactive IPython session.

For more details, see the `Dask distributed IPython docs
<https://distributed.dask.org/en/latest/ipython.html>`_.
Manage environments
===================

It is critical that each of your dask workers uses the same set of
python packages and modules when executing your code, so that Dask
can function. Upon connecting with a distributed ``Client``, Dask
will automatically check the versions of some critical packages
(including Dask itself) and warn you of any mismatch.

Most functions you will run on Dask will require imports. This
is even true for passing any object which is not a python builtin -
the pickle serialisation method will save references to imported modules
rather than trying to send all of your source code.

You therefore must ensure that workers have access to all of the modules
you will need, and ideally with exactly the same versions.

Single-machine schedulers
`````````````````````````

If you are using the threaded scheduler, then you do not need to do
anything, since the workers are in the same process, and objects are
simply shared rather than serialised and deserialised.

Similarly, if you use the multiprocessing scheduler, new processes
will be copied from, or launched in the same way as the original process,
so you only need make sure that you have not changed environment variables
related to starting up python and importing
code (such as PATH, PYTHONPATH, ``sys.path``).

If you are using the distributed scheduler on a single machine, this is roughly
equivalent to using the multiprocessing scheduler, above, if you are launching
with ``Client(...)`` or ``LocalCluster(...)``.

However, if you are launching your workers from the command line, then you must
ensure that you are running in the same environment (virtualenv, pipenv or conda).

The rest of this page concerns only distributed clusters.

Maintain consistent environments
````````````````````````````````

If you manage your environments yourself, then setting up module consistency
can be as simple as creating environments from the same pip or conda specification
on each machine. You should consult the documentation for ``pip``, ``pipenv``
and ``conda``, whichever you normally use. You will normally want to be as specific
about package versions as possible, and distribute the same environment file to
workers before installation.

However, other common ways to distribute an environment directly, rather than build it
in-place, include:

- docker images, where the environment has been built into the image; this is the
  normal route when you are running on infrastructure enabled by docker, such as
  kubernetes
- `conda-pack`_ is a tool for bundling existing conda environments, so they can be
  relocated to other machines. This tool was specifically created for dask on YARN/hadoop
  clusters, but could be used elsewhere
- shared filesystem, e.g., NFS, that can be seen by all machines. Note that importing
  python modules is fairly IO intensive, so your server needs to be able to handle
  many requests
- cluster install method (e.g., `parcels`_): depending on your infrastructure, there may be
  ways to install specific binaries to all workers in a cluster.

.. _conda-pack: https://conda.github.io/conda-pack/
.. _parcels: https://docs.cloudera.com/documentation/enterprise/latest/topics/cm_ig_parcels.html

Temporary installations
```````````````````````
The worker plugin ``distributed.diagnostics.plugin.PipInstall`` allows you to
run pip installation commands on your workers, and optionally have them restart
upon success. Please read the plugin documentation to see how to use this.

Objects in ``__main__``
```````````````````````

Objects that you create without any reference to modules, such as classes that
are defined right in the repl or notebook cells, are not pickled with reference to
any imported module. You can redefine such objects, and Dask will serialise them
completely, including the source code. This can be a good way to try new things
while working with a distributed setup.

Send Source
```````````

Particularly during development, you may want to send files directly to workers
that are already running.

You should use ``client.upload_file`` in these cases.
For more detail, see the `API docs`_ and a
StackOverflow question
`"Can I use functions imported from .py files in Dask/Distributed?"`__
This function supports both standalone file and setuptools's ``.egg`` files
for larger modules.

__ http://stackoverflow.com/questions/39295200/can-i-use-functions-imported-from-py-files-in-dask-distributed
.. _API docs: https://distributed.readthedocs.io/en/latest/api.html#distributed.executor.Executor.upload_file
Setup adaptive deployments
==========================

Motivation
----------

Most Dask deployments are static with a single scheduler and a fixed number of
workers.  This results in predictable behavior, but is wasteful of resources in
two situations:

1.  The user may not be using the cluster, or perhaps they are busy
    interpreting a recent result or plot, and so the workers sit idly,
    taking up valuable shared resources from other potential users
2.  The user may be very active, and is limited by their original allocation.

Particularly efficient users may learn to manually add and remove workers
during their session, but this is rare.  Instead, we would like the size of a
Dask cluster to match the computational needs at any given time.  This is the
goal of the *adaptive deployments* discussed in this document.

|

.. image:: ../images/dask-adaptive.svg
   :alt: Dask adaptive scaling
   :align: center
   :scale: 40%

|

These are particularly helpful for interactive workloads, which are characterized by long
periods of inactivity interrupted with short bursts of heavy activity.
Adaptive deployments can result in both faster analyses that give users much
more power, but with much less pressure on computational resources.

.. raw:: html

   <iframe width="560"
           height="315"
           src="https://www.youtube.com/embed/dViyEqOMA8U"
           style="margin: 0 auto 20px auto; display: block;"
           frameborder="0"
           allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
           allowfullscreen></iframe>


Adaptive
--------

To make setting up adaptive deployments easy, some Dask deployment solutions
offer an ``.adapt()`` method.  Here is an example with
`dask_kubernetes.KubeCluster
<https://kubernetes.dask.org/en/latest/kubecluster.html>`_.

.. code-block:: python

   from dask_kubernetes import KubeCluster

   cluster = KubeCluster()
   cluster.adapt(minimum=0, maximum=100)  # scale between 0 and 100 workers

For more keyword options, see the Adaptive class below:

.. currentmodule:: distributed.deploy

.. autosummary::
   Adaptive


Dependence on a Resource Manager
--------------------------------

The Dask scheduler does not know how to launch workers on its own. Instead, it
relies on an external resource scheduler like Kubernetes above, or
Yarn, SGE, SLURM, Mesos, or some other in-house system (see :doc:`how to deploy Dask
clusters <../how-to/deploy-dask-clusters>` for options).  In order to use adaptive deployments, you
must provide some mechanism for the scheduler to launch new workers.  Typically,
this is done by using one of the solutions listed in the :doc:`how to deploy Dask
clusters <../how-to/deploy-dask-clusters>`, or by subclassing from the Cluster superclass and
implementing that API.

.. autosummary::
   Cluster


Scaling Heuristics
------------------

The Dask scheduler tracks a variety of information that is useful to correctly
allocate the number of workers:

1.  The historical runtime of every function and task that it has seen,
    and all of the functions that it is currently able to run for users
2.  The amount of memory used and available on each worker
3.  Which workers are idle or saturated for various reasons, like the presence
    of specialized hardware

From these, it is able to determine a target number of workers by dividing the
cumulative expected runtime of all pending tasks by the ``target_duration``
parameter (defaults to five seconds).  This number of workers serves as a
baseline request for the resource manager.  This number can be altered for a
variety of reasons:

1.  If the cluster needs more memory, then it will choose either the target
    number of workers or twice the current number of workers (whichever is
    larger)
2.  If the target is outside of the range of the minimum and maximum values,
    then it is clipped to fit within that range

Additionally, when scaling down, Dask preferentially chooses those workers that
are idle and have the least data in memory.  It moves that data to other
machines before retiring the worker.  To avoid rapid cycling of the cluster up
and down in size, we only retire a worker after a few cycles have gone by where
it has consistently been a good idea to retire it (controlled by the
``wait_count`` and ``interval`` parameters).


API
---

.. autoclass:: Adaptive
.. autoclass:: Cluster
Setup Prometheus monitoring
===========================

Prometheus_ is a widely popular tool for monitoring and alerting a wide variety of systems. Dask.distributed exposes
scheduler and worker metrics in a prometheus text based format. Metrics are available at ``http://scheduler-address:8787/metrics`` when the ``prometheus_client`` package has been installed.

.. _Prometheus: https://prometheus.io

Available metrics are as following

+---------------------------------------------+------------------------------------------------+-----------+--------+
| Metric name                                 | Description                                    | Scheduler | Worker |
+=========================+===================+================================================+===========+========+
| python_gc_objects_collected_total           | Objects collected during gc.                   |    Yes    |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| python_gc_objects_uncollectable_total       | Uncollectable object found during GC.          |    Yes    |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| python_gc_collections_total                 | Number of times this generation was collected. |    Yes    |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| python_info                                 | Python platform information.                   |    Yes    |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_scheduler_workers                      | Number of workers connected.                   |    Yes    |        |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_scheduler_clients                      | Number of clients connected.                   |    Yes    |        |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_scheduler_tasks                        | Number of tasks at scheduler.                  |    Yes    |        |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_worker_tasks                           | Number of tasks at worker.                     |           |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_worker_connections                     | Number of task connections to other workers.   |           |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_worker_threads                         | Number of worker threads.                      |           |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_worker_latency_seconds                 | Latency of worker connection.                  |           |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_worker_tick_duration_median_seconds    | Median tick duration at worker.                |           |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_worker_task_duration_median_seconds    | Median task runtime at worker.                 |           |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+
| dask_worker_transfer_bandwidth_median_bytes | Bandwidth for transfer at worker in Bytes.     |           |  Yes   |
+---------------------------------------------+------------------------------------------------+-----------+--------+

How To...
=========

This section contains snippets and suggestions about how to perform different actions
using Dask. If you have an idea of a how-to that we should add, please
`make a suggestion <https://github.com/dask/dask/tree/main/docs/source/how-to>`_!

.. Articles in this section should be short and not contain much explanation.

.. toctree::
   :caption: How To...
   :maxdepth: 1
   :glob:

   *
   Use GPUs <../gpu.rst>
Connect to remote data
======================

Dask can read data from a variety of data stores including local file systems,
network file systems, cloud object stores, and Hadoop.  Typically this is done
by prepending a protocol like ``"s3://"`` to paths used in common data access
functions like ``dd.read_csv``:

.. code-block:: python

   import dask.dataframe as dd
   df = dd.read_csv('s3://bucket/path/to/data-*.csv')
   df = dd.read_parquet('gcs://bucket/path/to/data-*.parq')

   import dask.bag as db
   b = db.read_text('hdfs://path/to/*.json').map(json.loads)

Dask uses `fsspec`_ for local, cluster and remote data IO. Other file interaction, such
as loading of configuration, is done using ordinary python method.

The following remote services are well supported and tested against the main
codebase:

- **Local or Network File System**: ``file://`` - the local file system, default in the absence of any protocol.

- **Hadoop File System**: ``hdfs://`` - Hadoop Distributed File System, for resilient, replicated
  files within a cluster. This uses PyArrow_ as the backend.

- **Amazon S3**: ``s3://`` - Amazon S3 remote binary store, often used with Amazon EC2,
  using the library s3fs_.

- **Google Cloud Storage**: ``gcs://`` or ``gs://`` - Google Cloud Storage, typically used with Google Compute
  resource using gcsfs_.

- **Microsoft Azure Storage**: ``adl://``, ``abfs://`` or ``az://`` - Microsoft Azure Storage using adlfs_.

- **HTTP(s)**: ``http://`` or ``https://`` for reading data directly from HTTP web servers.

`fsspec`_ also provides other file systems that may be of interest to Dask users, such as
ssh, ftp, webhdfs and dropbox. See the documentation for more information.

When specifying a storage location, a URL should be provided using the general
form ``protocol://path/to/data``.  If no protocol is provided, the local
file system is assumed (same as ``file://``).

.. _fsspec: https://filesystem-spec.readthedocs.io/
.. _s3fs: https://s3fs.readthedocs.io/
.. _adlfs: https://github.com/dask/adlfs
.. _gcsfs: https://gcsfs.readthedocs.io/en/latest/
.. _PyArrow: https://arrow.apache.org/docs/python/

Lower-level details on how Dask handles remote data is described
below in the Internals section

Optional Parameters
-------------------

Two methods exist for passing parameters to the backend file system driver:
extending the URL to include username, password, server, port, etc.; and
providing ``storage_options``, a dictionary of parameters to pass on. The
second form is more general, as any number of file system-specific options
can be passed.

Examples:

.. code-block:: python

   df = dd.read_csv('hdfs://user@server:port/path/*.csv')

   df = dd.read_parquet('s3://bucket/path',
                        storage_options={'anon': True, 'use_ssl': False})

Details on how to provide configuration for the main back-ends
are listed next, but further details can be found in the documentation pages of the
relevant back-end.

Each back-end has additional installation requirements and may not be available
at runtime. The dictionary ``fsspec.registry`` contains the
currently imported file systems. To see which backends ``fsspec`` knows how
to import, you can do

.. code-block:: python

    from fsspec.registry import known_implementations
    known_implementations

Note that some backends appear twice, if they can be referenced with multiple
protocol strings, like "http" and "https".

Local File System
-----------------

Local files are always accessible, and all parameters passed as part of the URL
(beyond the path itself) or with the ``storage_options``
dictionary will be ignored.

This is the default back-end, and the one used if no protocol is passed at all.

We assume here that each worker has access to the same file system - either
the workers are co-located on the same machine, or a network file system
is mounted and referenced at the same path location for every worker node.

Locations specified relative to the current working directory will, in
general, be respected (as they would be with the built-in python ``open``),
but this may fail in the case that the client and worker processes do not
necessarily have the same working directory.

Hadoop File System
------------------

The Hadoop File System (HDFS) is a widely deployed, distributed, data-local file
system written in Java. This file system backs many clusters running Hadoop and
Spark. HDFS support can be provided by PyArrow_.

By default, the back-end attempts to read the default server and port from
local Hadoop configuration files on each node, so it may be that no
configuration is required. However, the server, port, and user can be passed as
part of the url: ``hdfs://user:pass@server:port/path/to/data``, or using the
``storage_options=`` kwarg.


Extra Configuration for PyArrow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following additional options may be passed to the ``PyArrow`` driver via
``storage_options``:

    - ``host``, ``port``, ``user``: Basic authentication
    - ``kerb_ticket``: Path to kerberos ticket cache

PyArrow's ``libhdfs`` driver can also be affected by a few environment
variables. For more information on these, see the `PyArrow documentation`_.

.. _PyArrow documentation: https://arrow.apache.org/docs/python/filesystems_deprecated.html#hadoop-file-system-hdfs


Amazon S3
---------

Amazon S3 (Simple Storage Service) is a web service offered by Amazon Web
Services.

The S3 back-end available to Dask is s3fs_, and is importable when Dask is
imported.

Authentication for S3 is provided by the underlying library boto3. As described
in the `auth docs`_, this could be achieved by placing credentials files in one
of several locations on each node: ``~/.aws/credentials``, ``~/.aws/config``,
``/etc/boto.cfg``, and ``~/.boto``. Alternatively, for nodes located
within Amazon EC2, IAM roles can be set up for each node, and then no further
configuration is required. The final authentication option for user
credentials can be passed directly in the URL
(``s3://keyID:keySecret/bucket/key/name``) or using ``storage_options``. In
this case, however, the key/secret will be passed to all workers in-the-clear,
so this method is only recommended on well-secured networks.

.. _auth docs: https://boto3.amazonaws.com/v1/documentation/api/latest/guide/configuration.html

The following parameters may be passed to s3fs using ``storage_options``:

    - anon: Whether access should be anonymous (default False)

    - key, secret: For user authentication

    - token: If authentication has been done with some other S3 client

    - use_ssl: Whether connections are encrypted and secure (default True)

    - client_kwargs: Dict passed to the `boto3 client`_, with keys such
      as `region_name` or `endpoint_url`. Notice: do not pass the `config`
      option here, please pass it's content to `config_kwargs` instead.

    - config_kwargs: Dict passed to the `s3fs.S3FileSystem`_, which passes it to
      the `boto3 client's config`_ option.

    - requester_pays: Set True if the authenticated user will assume transfer
      costs, which is required by some providers of bulk data

    - default_block_size, default_fill_cache: These are not of particular
      interest to Dask users, as they concern the behaviour of the buffer
      between successive reads

    - kwargs: Other parameters are passed to the `boto3 Session`_ object,
      such as `profile_name`, to pick one of the authentication sections from
      the configuration files referred to above (see `here`_)

.. _boto3 client: https://boto3.amazonaws.com/v1/documentation/api/latest/reference/core/session.html#boto3.session.Session.client
.. _boto3 Session: https://boto3.amazonaws.com/v1/documentation/api/latest/reference/core/session.html
.. _here: https://boto3.amazonaws.com/v1/documentation/api/latest/guide/credentials.html#shared-credentials-file
.. _s3fs.S3FileSystem: https://s3fs.readthedocs.io/en/latest/api.html#s3fs.core.S3FileSystem
.. _boto3 client's config: https://botocore.amazonaws.com/v1/documentation/api/latest/reference/config.html

Using Other S3-Compatible Services
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By using the `endpoint_url` option, you may use other s3-compatible services,
for example, using `AlibabaCloud OSS`:

.. code-block:: python

    dask_function(...,
        storage_options={
            "key": ...,
            "secret": ...,
            "client_kwargs": {
                "endpoint_url": "http://some-region.some-s3-compatible.com",
            },
            # this dict goes to boto3 client's `config`
            #   `addressing_style` is required by AlibabaCloud, other services may not
            "config_kwargs": {"s3": {"addressing_style": "virtual"}},
        })

Google Cloud Storage
--------------------

Google Cloud Storage is a RESTful online file storage web service for storing
and accessing data on Google's infrastructure.

The GCS back-end is identified by the
protocol identifiers ``gcs`` and ``gs``, which are identical in their effect.

Multiple modes of authentication are supported. These options should be
included in the ``storage_options`` dictionary as ``{'token': ..}`` submitted with your call
to a storage-based Dask function/method. See the `gcsfs`_ documentation for further
details.


General recommendations for distributed clusters, in order:

- use ``anon`` for public data
- use ``cloud`` if this is available
- use `gcloud`_ to generate a JSON file, and distribute this to all workers, and
  supply the path to the file

- use gcsfs directly with the ``browser`` method to generate a token cache file
  (``~/.gcs_tokens``) and distribute this to all workers, thereafter using method ``cache``

.. _gcloud: https://cloud.google.com/sdk/docs/

A final suggestion is shown below, which may be the fastest and simplest for authenticated access (as
opposed to anonymous), since it will not require re-authentication. However, this method
is not secure since credentials will be passed directly around the cluster. This is fine if
you are certain that the cluster is itself secured. You need to create a ``GCSFileSystem`` object
using any method that works for you and then pass its credentials directly:

.. code-block:: python

    gcs = GCSFileSystem(...)
    dask_function(..., storage_options={'token': gcs.session.credentials})

Microsoft Azure Storage
-----------------------

Microsoft Azure Storage is comprised of Data Lake Storage (Gen1) and Blob Storage (Gen2).
These are identified by the protocol identifiers ``adl`` and ``abfs``, respectively,
provided by the adlfs_ back-end.

Authentication for ``adl`` requires ``tenant_id``, ``client_id`` and ``client_secret``
in the ``storage_options`` dictionary.

Authentication for ``abfs`` requires ``account_name`` and ``account_key`` in ``storage_options``.

HTTP(S)
-------

Direct file-like access to arbitrary URLs is available over HTTP and HTTPS. However,
there is no such thing as ``glob`` functionality over HTTP, so only explicit lists
of files can be used.

Server implementations differ in the information they provide - they may or may
not specify the size of a file via a HEAD request or at the start of a download -
and some servers may not respect byte range requests. The HTTPFileSystem therefore
offers best-effort behaviour: the download is streamed but, if more data is seen
than the configured block-size, an error will be raised. To be able to access such
data you must read the whole file in one shot (and it must fit in memory).

Using a block size of 0 will return normal ``requests`` streaming file-like objects,
which are stable, but provide no random access.

Developer API
-------------

The prototype for any file system back-end can be found in
``fsspec.spec.AbstractFileSystem``. Any new implementation should provide the
same API, or directly subclass, and make itself available as a protocol to Dask. For example, the
following would register the protocol "myproto", described by the implementation
class ``MyProtoFileSystem``. URLs of the form ``myproto://`` would thereafter
be dispatched to the methods of this class:

.. code-block:: python

   fsspec.registry['myproto'] = MyProtoFileSystem

However, it would be better to submit a PR to ``fsspec`` to include the class in
the ``known_implementations``.

Internals
---------

Dask contains internal tools for extensible data ingestion in the
``dask.bytes`` package and using `fsspec`_.
.  *These functions are developer-focused rather than for
direct consumption by users.  These functions power user facing functions like*
``dd.read_csv`` *and* ``db.read_text`` *which are probably more useful for most
users.*


.. currentmodule:: dask.bytes

.. autosummary::
   read_bytes
   open_files

These functions are extensible in their output formats (bytes, file objects),
their input locations (file system, S3, HDFS), line delimiters, and compression
formats.

Both functions are *lazy*, returning either
pointers to blocks of bytes (``read_bytes``) or open file objects
(``open_files``).  They can handle different storage backends by prepending
protocols like ``s3://`` or ``hdfs://`` (see below). They handle compression formats
listed in ``fsspec.compression``, some of which may require additional packages
to be installed.

These functions are not used for all data sources.  Some data sources like HDF5
are quite particular and receive custom treatment.

Delimiters
^^^^^^^^^^

The ``read_bytes`` function takes a path (or globstring of paths) and produces
a sample of the first file and a list of delayed objects for each of the other
files.  If passed a delimiter such as ``delimiter=b'\n'``, it will ensure that
the blocks of bytes start directly after a delimiter and end directly before a
delimiter.  This allows other functions, like ``pd.read_csv``, to operate on
these delayed values with expected behavior.

These delimiters are useful both for typical line-based formats (log files,
CSV, JSON) as well as other delimited formats like Avro, which may separate
logical chunks by a complex sentinel string. Note that the delimiter finding
algorithm is simple, and will not account for characters that are escaped,
part of a UTF-8 code sequence or within the quote marks of a string.

Compression
^^^^^^^^^^^

These functions support widely available compression technologies like ``gzip``,
``bz2``, ``xz``, ``snappy``, and ``lz4``.  More compressions can be easily
added by inserting functions into dictionaries available in the
``fsspec.compression`` module.  This can be done at runtime and need not be
added directly to the codebase.

However, most compression technologies like ``gzip`` do not support efficient
random access, and so are useful for streaming ``open_files`` but not useful for
``read_bytes`` which splits files at various points.

API
^^^

.. autofunction:: read_bytes
.. autofunction:: open_files
