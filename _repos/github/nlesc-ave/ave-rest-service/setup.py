from setuptools import setup, find_packages

exec(open('avedata/version.py').read())

setup(
    name='avedata',
    packages=find_packages(),
    version=__version__,
    include_package_data=True,
    install_requires=[
        'flask',
        'click',
        'connexion',
        'whoosh',
        'pybedtools',
        'cyvcf2',
        'numpy',
        'scipy',
    ],
    package_data={
      'avedata': ['schema.sql', 'swagger.yml']
    },
    entry_points='''
        [console_scripts]
        avedata=avedata.commands:cli
    ''',
    license='Apache',
    classifiers=[
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Framework :: Flask',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Environment :: Web Environment',
    ],
)
