from setuptools import setup, find_packages
 
setup(name='komics',
      version='1.1.8',
      description='Kinetoplast genOMICS',
      long_description='Assembly and circularization of mitochondrial minicircles',
      author='Frederik Van den Broeck',
      author_email='fvandenbroeck@kuleuven.be',
      url='http://github.com/frebio/komics',
      keywords='megahit blast assembly circularization',
      license='GPLv3',
      scripts=['bin/komics'],
      packages=find_packages(exclude=['docs', 'tests*', 'data*']),
      package_data={'komics': ['data/*']},
      include_package_data=True,
      install_requires=[
	'biopython',
      ],
      #dependency_links=['https://github.com/sanger-pathogens/Fastaq.git'],
      classifiers=[
      	'Development Status :: 2 - Pre-Alpha',
      	'Topic :: Scientific/Engineering :: Bio-Informatics',
      	'Programming Language :: Python :: 3.7',
      	'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
      zip_safe=False,
)
