from setuptools import setup

setup(
    name='rMLST-API',
    version='0.0.7',
    py_modules=['src.rmlst_api.*'],
    install_requires=[
        'Click',
    ],
    entry_points={
        'console_scripts': [
            'rmlst-api = src.rmlst_api.cli:run_all',
        ],
    },
)