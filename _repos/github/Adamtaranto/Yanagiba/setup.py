from setuptools import setup

pypi_classifiers = [
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    "Topic :: Software Development :: Libraries :: Python Modules",
    'License :: OSI Approved :: MIT License',
]

install_requires = [
	"nanomath>=0.13.0",
  "pandas>=0.20.3",
  'biopython>=1.70',
]

desc = """Filter short or low quality Oxford Nanopore reads which have been basecalled with Albacore."""

setup(name='yanagiba',
      version='1.0.0',
      description=desc,
      url='https://github.com/Adamtaranto/Yanagiba',
      author='Adam Taranto',
      author_email='adam.taranto@anu.edu.au',
      license='MIT',
      packages=['yanagiba'],
      classifiers=pypi_classifiers,
      keywords=["Albacore","Nanopore","basecalling","genome","DNA","sequencing"],
      install_requires=install_requires,
      include_package_data=True,
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'yanagiba=yanagiba.cmd_line:main',
        ],
    },
    )