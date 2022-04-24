"""
metapredict
A protein disorder predictor based on a BRNN (IDP-Parrot) trained on the consensus disorder values from 8 disorder predictors from 12 proteomes.
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
pytest_runner = ['pytest-runner'] if needs_pytest else []

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    # Self-descriptive entries which should always be present
    name='metapredict',
    author='Ryan Emenecker - Holehouse Lab - WUSM',
    author_email='remenecker@wustl.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url = 'https://github.com/idptools/metapredict.git',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is installed
    # Handled automatically by setuptools. Use 'exclude' to prevent some specific
    # subpackage(s) from being added, if needed
    packages=find_packages(),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,

    # Allows `setup.py test` to work correctly with pytest
    setup_requires=[] + pytest_runner,
    scripts=['scripts/metapredict-graph-disorder',
             'scripts/metapredict-predict-disorder',
             'scripts/metapredict-quick-graph',
             'scripts/metapredict-quick-predict',
             'scripts/metapredict-uniprot',
             'scripts/metapredict-predict-idrs',
             'scripts/metapredict-graph-pLDDT',
             'scripts/metapredict-predict-pLDDT',
             'scripts/metapredict-name'],             
         

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # url='http://www.my_package.com',  # Website
    install_requires=[
            'torch',
            'numpy',
            'matplotlib',
            'protfasta',
            'scipy',
            'urllib3',
            'alphaPredict==1.0'],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
