from distutils.core import setup
from setuptools import find_packages

# Setup optional dependencies
optional_dependencies = {
    'docs': [
        'jupyter',
        'jupyter-book',
        'nbconvert',
        'notebook',
        'ipython',
        'numpydoc>=0.9.1',
        'redbaron',
        'tabulate',
        'jaxlib',
        'jax'
    ],
    'notebooks': [
        'notebook',
        'tabulate',
        'ipython'
    ],
    'test': [
        'pep8',
        'testflo>=1.3.6',
        'matplotlib',
        'numpydoc>=1.1',
    ]
}

# Add an optional dependency that concatenates all others
optional_dependencies['all'] = sorted([
    dependency
    for dependencies in optional_dependencies.values()
    for dependency in dependencies
])


setup(name='dymos',
    version='1.4.1-dev',
    description='Open-Source Optimization of Dynamic Multidisciplinary Systems',
    long_description='''
Dymos is a framework for the simulation and optimization of dynamical systems within the OpenMDAO Multidisciplinary Analysis and Optimization environment.
Dymos leverages implicit and explicit simulation techniques to simulate generic dynamic systems of arbitary complexity.

The software has two primary objectives:
-   Provide a generic ODE integration interface that allows for the analysis of dynamical systems.
-   Allow the user to solve optimal control problems involving dynamical multidisciplinary systems.''',
    long_description_content_type='text/markdown',
    url='https://github.com/OpenMDAO/dymos',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    license='Apache License',
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        'openmdao>=3.17.0',
        'numpy>=1.17',
        'scipy>=1.2'
    ],
    extras_require=optional_dependencies,
    zip_safe=False,
    package_data={'dymos.examples.aircraft_steady_flight.aero': ['data/CRM_aero_inputs.dat', 'data/CRM_aero_outputs.dat']},
)
