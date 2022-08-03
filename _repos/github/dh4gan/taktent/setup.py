from setuptools import setup, find_packages

setup(name='taktent',
      version='1.13',
      description='A framework for conducting agent-based simulations of SETI',
      long_description=open('README.md').read(),
      url='http://github.com/dh4gan/taktent',
      author='Duncan Forgan',
      author_email='dh4gan@gmail.com',
      license='GPL-3.0',
      install_requires=['numpy','matplotlib'],
      packages=['taktent','taktent.agents', 'taktent.strategies','taktent.populations'],
      zip_safe=False)
