from setuptools import setup

setup(name='PyEscape',
      version='2.0',
      description='Library used for simulating narrow escape problems',
      url='https://github.com/AoifeHughes/narrow_escape',
      author='Aoife Hughes',
      author_email='Aoife.hughes@jic.ac.uk',
      license='MIT',
      packages=['PyEscape'],
      install_requires=['numpy',
                        'matplotlib',
                        'scipy',
                        'plotly',
                        'tqdm'],
      entry_points={
          'console_scripts': [
              'PyEscape = PyEscape.__main__:main'
          ]
      },
      zip_safe=True)
