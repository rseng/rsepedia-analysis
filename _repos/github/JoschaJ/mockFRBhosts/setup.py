from setuptools import setup


setup(name='mockFRBhosts',
      version='1.0.0',
      description='Package to simulate FRB hosts and their follow-up.',
      author='Joscha N. Jahns',
      author_email='jjahns@mpifr-bonn.mpg.de',
      url='https://github.com/JoschaJ/mockFRBhosts',
      packages=['mockFRBhosts'],
      license='MIT',
      python_requires='>=3.6',  # because of f-strings
      install_requires=[
          'numpy<=1.22',  # because of pymc3 and frbpoppy using np.warnings
          'scipy',
          'pandas',
          'matplotlib',
          'seaborn',
          'astropy',
          'corner',
          ]
      )
