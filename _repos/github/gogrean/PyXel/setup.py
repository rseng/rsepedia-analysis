from setuptools import setup

setup(name='pyxel',
      version='0.1.0',
      description='Astronomical X-ray Image Analysis Package',
      author='Georgiana Ogrean',
      author_email='georgiana.ogrean@gmail.com',
      url='https://github.com/gogrean/PyXel',
      packages=['pyxel'],
      install_requires=[
          'astropy>=1.1.2',
          'matplotlib>=1.5.1',
          'numpy>=1.11',
          'scipy>=0.17',
          'emcee>=2.1',
          'corner>=1.0.2',
          'tabulate>=0.7.5'
      ]
     )
