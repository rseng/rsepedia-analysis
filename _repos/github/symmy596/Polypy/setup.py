__author__ = "Adam R. Symington"
__copyright__ = "Copyright Adam R.Symington (2019)"
__version__ = "1.0.0"
__maintainer__ = "Adam R. Symington"
__email__ = "ars44s@bath.ac.uk"
__date__ = "27/02/2021"

from setuptools import setup
import os

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='polypy',
        version='1.0.0',
        description='Molecular Dynamics analysis',
        long_description=open(os.path.join(module_dir, 'README.md')).read(),
        long_description_content_type='text/markdown', 
        url='https://github.com/symmy596/Polypy',
        author='Adam R. Symington',
        author_email='ars44@bath.ac.uk',
        license='MIT license',
        packages=['polypy'],
        zip_safe=False,
        install_requires=['scipy','numpy'],
        classifiers=['Programming Language :: Python',
                     'Development Status :: 5 - Production/Stable',
                     'Intended Audience :: Science/Research',
                     'Operating System :: OS Independent',
                     'Topic :: Scientific/Engineering']
    )
