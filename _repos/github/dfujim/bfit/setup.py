
# install numpy and cython dependencies needed for the rest of the setup script
from setuptools import dist
dist.Distribution().fetch_build_eggs(['cython>=0.28', 'numpy>=1.19'])

import setuptools
from distutils.core import Extension
from Cython.Build import cythonize
import numpy, os, sys
from os.path import join

with open("README.md", "r", encoding = "utf8") as fh:
    long_description  =  fh.read()

# get setup variables
variables  =  {}
with open(os.path.join('bfit', 'global_variables.py')) as fid:
    exec(fid.read(), variables)
    
__version__  =  variables['__version__']
__src__  =  variables['__src__']

# get libraries
libraries = []
if sys.platform in ('unix', 'darwin'):
    libraries.append('m')

# module extension
ext  =  Extension("bfit.fitting.integrator",
                sources = [join(__src__, "integrator.pyx"),
                           join(__src__, "FastNumericalIntegration_src", "integration_fns.cpp")],
                language = "c++",             # generate C++ code                        
                include_dirs = [join(__src__, "FastNumericalIntegration_src"), 
                                numpy.get_include()],
                libraries = libraries,
                extra_compile_args = ["-ffast-math", '-O3']
                )

setuptools.setup(
    name = "bfit",
    version = __version__,
    author = "Derek Fujimoto",
    author_email = "fujimoto@phas.ubc.ca",
    description = "β-NMR and β-NQR Data Analysis",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/dfujim/bfit",
    packages = setuptools.find_packages(),
    classifiers = [
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",     
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Cython",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Development Status :: 5 - Production/Stable",
    ],
    install_requires = ['cython >= 0.28', 'numpy >= 1.19', 'tqdm >= 4.25.0',
                      'bdata >= 6.8.2', 'matplotlib >= 2.2.4', 'pandas >= 0.23.0',
                      'pyyaml >= 5.1', 'scipy >= 1.2.0', 'iminuit >= 2.6.1', 
                      'requests >= 2.25.0', 'argparse >= 1.4.0', 'pytest >= 4.5.0',
                      'wheel>=0.34', 'jax>=0.2.17', 'jaxlib>=0.1.69'],
    package_data = {'': ['data']},
    entry_points = {'console_scripts':['bfit = bfit:main']},
    include_package_data = True,
    ext_modules  =  cythonize([ext], include_path = [numpy.get_include()]),
)
