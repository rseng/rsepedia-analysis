from setuptools import setup, find_packages

setup(
    name="gims",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "flask==1.1.2",
        "ase==3.21.1",
        "spglib==1.15.0",
        "sphinx",
        "sphinx_rtd_theme",
        "pytest-cov==2.8.1",
        "pytest-remotedata==0.3.2",
    ],
    extras_require={
        'dev': [
            "sphinx",
            "sphinx_rtd_theme",
            "pytest-cov==2.8.1",
            "pytest-remotedata==0.3.2",
        ]
    },
    version="1.0.9",
)
