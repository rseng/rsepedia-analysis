# setup.py

from codecs import open
from os import path
import sys
from setuptools import setup

here = path.abspath(path.dirname(__file__))

with open(path.join(here, '_version.py')) as version_file:
    exec(version_file.read())

with open(path.join(here, 'README.md')) as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'CHANGELOG.md')) as changelog_file:
    changelog = changelog_file.read()

desc = readme + '\n\n' + changelog
try:
    import pypandoc
    long_description = pypandoc.convert_text(desc, 'rst', format='md')
    with open(path.join(here, 'README.rst'), 'w') as rst_readme:
        rst_readme.write(long_description)
except (ImportError, OSError, IOError):
    long_description = desc

install_requires = [
    'pandas',
    'apyori',
]

tests_require = [
    'pytest',
    'pytest-cov',
]

needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
setup_requires = ['pytest-runner'] if needs_pytest else []

setup(
    name='autofunc',
    version=__version__,
    description='Automatic Functional Representation Generation',
    long_description=long_description,
    author='Alex Mikes',
    author_email='alexmikes@gmail.com',
    url='https://github.com/DesignEngrLab',
    classifiers=[
        'License :: OSI Approved :: BSD License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
    ],
    license='BSD-3-Clause',
    install_requires=install_requires,
    tests_require=tests_require,
    python_requires='>=3',
    setup_requires=setup_requires,
    zip_safe=False,
    packages=['autofunc', 'autofunc.tests'],
    package_dir={
        'autofunc': 'autofunc',
        'autofunc.tests': 'autofunc/tests',
        },
    include_package_data=True,
)