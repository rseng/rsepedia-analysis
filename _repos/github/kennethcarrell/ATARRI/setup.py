from setuptools import find_packages, setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='ATARRI',
    version='v1.1',
    description='A TESS Archive RR Lyrae Classifier',
    python_requires='==3.8.*',
    long_description=readme(),
    platforms=['any'],
    packages=['ATARRI'],
    url='https://github.com/kennethcarrell/ATARRI',
    license='GPL-3.0',
    author='Kenneth Carrell',
    author_email='kennethcarrell@gmail.com',
    install_requires=['numpy>=1.19', 'astropy>=4.2', 'lightkurve>=2.0', 'astroquery>=0.4', 'matplotlib>=3.3', 'scipy>=1.6'],
    keywords=['astronomy', 'periodic variables', 'light curves', 'RR Lyrae', 'TESS'],
    classifiers=['Intended Audience :: Science/Research', 'Topic :: Scientific/Engineering :: Astronomy', 'Development Status :: 5 - Production/Stable', 'Operating System :: OS Independent', 'Programming Language :: Python'],
)

