import setuptools

with open('README.md', 'r') as file:
    long_description = file.read()

setuptools.setup(
    name = 'richvalues',
    version = '1.2.0',
    license = 'BSD-3-Clause',
    author = 'Andrés Megías Toledano',
    description = 'Python library for dealing with uncertainties and upper/lower limits',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    packages = setuptools.find_packages('./'),
    url = 'https://github.com/andresmegias/richvalues',
    install_requires = ['numpy', 'pandas', 'scipy', 'matplotlib'],
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent'
    ]
)
