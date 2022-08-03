from setuptools import setup
from os import path


try:  # Python 2
    execfile(path.join("EDIunplugged", 'version.py'))
except:  # Python 3
    exec(open(path.join("EDIunplugged", 'version.py')).read())


# If Python3: Add "README.md" to setup. 
# Useful for PyPI (pip install transitleastsquares). Irrelevant for users using Python2
# try:
#
#     this_directory = path.abspath(path.dirname(__file__))
#     with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
#         long_description = f.read()
# except:
#     long_description = ' '

    
setup(name='EDIunplugged',
    version=EDI_VERSIONING,
    description='An algorithm to identify false positive transit signals using TLS information',
    url='https://github.com/jonzink/EDI_Vetter_unplugged',
    author='Jon Zink',
    author_email='jzink@astro.ucla.edu',
    license='MIT',
    packages=['EDIunplugged'],
    include_package_data=True,
    # package_data={'': ['*.csv', '*.cfg']},
    # entry_points = {'console_scripts': ['EDIunplugged=EDIunplugged.command_line:main'],},
    install_requires=[
        'astropy',
        'numpy',
        'scipy',
        'transitleastsquares',
        ]
)
