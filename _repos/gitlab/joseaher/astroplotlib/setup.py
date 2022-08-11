from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='astroplotlib',
      version='0.2.4',
      description='Python scripts to handle astronomical images',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='J. A. Hernandez-Jimenez',
      author_email='joseaher@gmail.com',
      url = "https://gitlab.com/joseaher/astrocubelib",
      packages=['plot_functions', 'plot_slit', 'plot_ellipse', 'cata', 'sky'],
      scripts = ['scripts/guicata.py', 'scripts/calculate_sky.py',
                 'scripts/mag.sh', 'scripts/guiregions.py', 'scripts/slit.py',
                 'scripts/calibration.py', 'scripts/guiellipse.py',
                 'scripts/elme_images.py']
      )
