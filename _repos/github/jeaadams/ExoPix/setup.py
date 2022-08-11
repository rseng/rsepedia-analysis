from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='ExoPix',
    version='0.1.0',
    description='ExoPix: Applying PSF Subtraction with pyKLIP to JWST',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/jeaadams/ExoPix.git',
    author='Jea Adams, Jason Wang',
    author_email='jadams21@amherst.edu',
    include_package_data = True,
    packages=find_packages(),
    classifiers=[
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.8',
        ],
    keywords='KLIP PSF Subtraction Exoplanets Astronomy JWST',
    install_requires=['numpy', 'scipy', 'pyklip', 'astropy', 'matplotlib', 'emcee', 'corner', 'nbsphinx']
    )
