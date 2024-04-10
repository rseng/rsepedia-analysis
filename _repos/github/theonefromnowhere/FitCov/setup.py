import setuptools
import os
import sys

long_description = 'FitCov'


package_basedir = os.path.abspath(os.path.dirname(__file__))
package_basename = 'FitCov'

sys.path.insert(0, os.path.join(package_basedir, package_basename))



if __name__ == '__main__':
    setuptools.setup(
        name="FitCov",
        version="1.0",
        author="Svyatoslav Trusov",
        author_email="strusov@lpnhe.in2p3.fr",
        description="Codes for fitted jackknife covariance",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/theonefromnowhere/ClusteringLibs",
        packages=[package_basename],
                #setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        python_requires='>=3.6',
        install_requires=['numpy','pycorr','iminuit'],
    )