from setuptools import setup, find_packages

setup(
    name="SeqPanther",
    version="0.0.1",
    install_requires=[
        'biopython', 'click', 'numpy', 'pandas', 'pyfaidx', 'pysam',
        'matplotlib'
    ],
    packages=find_packages(include=["seqPanther", "seqPanther.*"]),
    entry_points={
        'console_scripts': ['seqpanther=seqPanther.seqPanther:run'],
    },
    url="https://github.com/krisp-kwazulu-natal/seqPatcher",
    license="GPLv3",
    author="Anmol Kiran; San James Emmanuel",
    author_email="anmol.kiran@gmail.com;sanemmanueljames@gmail.com",
    description="A set of sequence manipulation tools",
)
