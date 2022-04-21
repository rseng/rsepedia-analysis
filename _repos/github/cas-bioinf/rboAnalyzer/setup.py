import sys
from setuptools import setup

# hack to support wide range of configurations
if sys.version_info.minor == 4:
    # tested on ubuntu 14.04
    # pip3 and setuptools need to be updated prior to installation (see below)
    # pip3 install --user --upgrade "pip>18, <19.2" "setuptools>42, <44"
    subset = ['matplotlib>2.2, <3', 'biopython>=1.72, <1.74a0', 'pandas>=0.22, <0.24a0', 'numpy>=1.14, <1.16a0', 'MarkupSafe>=0.23.0, <2.0.0a0', 'kiwisolver>=1.0.1, <1.2']
elif sys.version_info.minor == 5:
    # tested on ubuntu 16.04
    # pip3 install --user "pip>20, <20.2.3" "setuptools>49, <51"
    subset = ['matplotlib>2.2, <3.1', 'biopython>=1.72, !=1.74, <1.77a0', 'pandas>=0.22, <0.25a0', 'numpy>=1.14, <1.19a0']
elif sys.version_info.minor == 6:
    # tested on ubuntu 18.04
    subset = ['matplotlib>2.2, <3.4a0', 'biopython==1.77', 'pandas>=0.22, <1', 'numpy>=1.14, <1.20a0']
elif sys.version_info.minor == 7:
    # tested on ubuntu 19.04
    subset = ['matplotlib>2.2, <3.4a0', 'biopython==1.77', 'pandas>=0.22, <1', 'numpy>=1.15, <1.21a0']
elif sys.version_info.minor == 8:
    # tested on ubuntu 20.04
    subset = ['matplotlib>3, <3.4a0', 'biopython==1.77', 'pandas>=0.22, <1', 'numpy>=1.18, <1.21a0']
else:
    subset = ['matplotlib', 'biopython>=1.72, !=1.74, <1.78a0', 'pandas>=0.22, <1', 'numpy>=1.14, <2']

with open('rna_blast_analyze/VERSION', 'r') as o:
    version = o.read().strip()

package_data = {
    'rna_blast_analyze': [
        'BR_core/config.txt',
        'BR_core/prediction_parameters.json',
        'BR_core/output/*',
        '3rd_party_source/RSEARCH_matrices/*',
        'docs/*',
        'VERSION'
    ]
}

setup(
    name='rboAnalyzer',
    version=version,
    description='Analyze BLAST output for RNA query',
    author='Marek Schwarz',
    author_email='marek.schwarz@biomed.cas.cz',
    entry_points={
        'console_scripts': [
            'rboAnalyzer = rna_blast_analyze.BA:main',
            'genomes_from_blast = rna_blast_analyze.download_blast_genomes:main'
        ],
    },
    install_requires=[
        'Jinja2>=2.9, <3',
    ] + subset,
    python_requires='>=3.4, <4',
    packages=['rna_blast_analyze', 'rna_blast_analyze.BR_core', 'rna_blast_analyze.BR_core.output'],
    package_data=package_data,
    include_package_data=True,
)
