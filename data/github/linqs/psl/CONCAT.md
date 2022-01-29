PSL
===

### Build Status
[![Master](https://travis-ci.org/linqs/psl.svg?branch=master)](https://travis-ci.org/linqs/psl)
[![Develop](https://travis-ci.org/linqs/psl.svg?branch=develop)](https://travis-ci.org/linqs/psl)

Probabilistic soft logic (PSL) is a machine learning framework for developing probabilistic models.
PSL models are easy to use and fast.
You can define models using a straightforward logical syntax and solve them with fast convex optimization.
PSL has produced state-of-the-art results in many areas spanning natural language processing, social-network analysis, knowledge graphs, recommender system, and computational biology.
More information about PSL is available at the [PSL homepage](https://psl.linqs.org).

Getting Started with PSL
------------------------

If you want to use PSL to build models, you probably do not need this source code.
Instead, visit the [Getting Started guide](https://psl.linqs.org/blog/2018/07/15/getting-started-with-psl.html) to learn how to create PSL projects that will automatically install a stable version of these libraries.

Installing PSL from Source
--------------------------

If you do want to install PSL from source, you can use [Maven](https://maven.apache.org/) 3.x.
In the top-level directory of the PSL source (which should be the same directory that holds this README), run:
```sh
mvn install
```

Citing PSL
----------

We hope you find PSL useful!
If you have, please consider citing PSL in any related publications as
```
@article{bach:jmlr17,
  Author = {Bach, Stephen H. and Broecheler, Matthias and Huang, Bert and Getoor, Lise},
  Journal = {Journal of Machine Learning Research (JMLR)},
  Title = {Hinge-Loss {M}arkov Random Fields and Probabilistic Soft Logic},
  Year = {2017}
}
```

Additional Resources
====================
- [PSL Homepage](https://psl.linqs.org)
- [PSL Examples](https://github.com/linqs/psl-examples)
- [API Reference](https://psl.linqs.org/api/)
- [PSL Source Repository](https://github.com/linqs/psl)
- [PSL Wiki](https://psl.linqs.org/wiki/)
- [Getting Started Guide](https://psl.linqs.org/blog/2018/07/15/getting-started-with-psl.html)
- [User Group](https://groups.google.com/forum/#!forum/psl-users)
## pslpython

A Python interface to the PSL SRL framework.

Instead of trying to fit a Python project into Maven conventions,
this package is generally formatted as a standard Python package.
However, special executions have been added to the following phases:
 - clean
    - `rm -r build dist pslpython.egg-info`
 - package
    - `python3 setup.py bdist_wheel`
 - install
    - `pip install --user --upgrade dist/pslpython-*.whl`
 - integration-test
    - `./run_tests.py`
 - deploy
    - `twine upload --repository-url https://test.pypi.org/legacy/ dist/pslpython-*.whl`
    - The `TWINE_USERNAME` and `TWINE_PASSWORD` environment variables MUST be set.

Instead of `compile` and `test`, `package` and `integration-test` are used.
This is because building the package relies on the jar from the `package` phase of psl-cli,
and the tests require the package to be built.

Also note that the star in the `install` phase means that the developer should clean before and install if versions have changed.

Because of the dependence on the parent's pom and psl-cli's jar, this package is only distributed as a wheel and not a source distribution.

Developers who are not using maven should make sure that the `psl-cli` package is built before building this package.
## pslpython

A Python interface to the PSL SRL framework.

PSL information can be found at: [psl.linqs.org](https://psl.linqs.org/).
Example usage can be in the PSL Examples repository:
 - [Stable](https://github.com/linqs/psl-examples)
 - [Development](https://github.com/linqs/psl-examples/tree/develop)
