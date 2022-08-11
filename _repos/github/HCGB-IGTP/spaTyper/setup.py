import setuptools
import glob

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

setuptools.setup(
    name="spaTyper",
    version="0.3.3",

    scripts=glob.glob('main/*'),
    author="Mitchell Sullivan; Jose F. Sanchez-Herrero",

    author_email="jfbioinformatics@gmail.com",
    description="Typing of Staphylococcus aureus protein A from assemblies",

    long_description_content_type="text/markdown",
    long_description=long_description_text,
    url="https://github.com/JFsanchezherrero/spa_typing",
    packages=setuptools.find_packages(),
    license='GPLv3',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ]
)
