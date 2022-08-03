from setuptools import setup


setup(
    name='yanosim',
    version='0.1',
    description=(
        'nanopore DRS read simulator'
    ),
    author='Matthew Parker',
    entry_points={
        'console_scripts': [
            'yanosim = yanosim.main:yanosim'
        ]
    },
    packages=[
        'yanosim',
    ],
    install_requires=[
        'numpy',
        'scipy',
        'click',
        'pysam',
    ],
)
