import os
from setuptools import setup, find_packages


directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(directory, 'version.py'), 'r') as f:
    exec(f.read())

setup(
    name='autostreamtree',
    version=__version__,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'autostreamtree=autostreamtree.cli:main',
            'networkDimensions=autostreamtree_scripts.networkDimensions:main',
            'streeToDendrogram=autostreamtree_scripts.streeToDendrogram:main',
        ],
    },
    install_requires=[
        "numpy",
        "pandas>=2.0",
        "networkx>=3.0",
        "seaborn",
        "matplotlib",
        "geopandas",
        "pyogrio",
        "momepy",
        "scipy",
        "scikit-learn",
        "mantel",
        "pysam",
        "sortedcontainers"
    ],
    extras_require={
        'dev': ['pytest']
    },
    package_data={
        'autostreamtree': ['data/*'],
    },
    include_package_data=True,
    author='Tyler Chafin',
    author_email='tyler.chafin@bioss.ac.uk',
    description='A package for fitting genetic distances to spatial networks',
    keywords='population genetics, genetic distance, river networks',
    url='https://github.com/tkchafin/autostreamtree',
)
