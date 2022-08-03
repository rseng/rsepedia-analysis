
from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='cosmic_kite',
    #url='https://github.com/jladan/package_demo',
    author='Martín de los Rios',
    author_email='martindelosrios13@gmail.com',
    # Needed to actually package something
    packages=['cosmic_kite'],
    package_dir={'cosmic_kite': 'src/'},
    # Needed for dependencies
    install_requires=['numpy','tensorflow','scikit-learn==0.22.2.post1'],
    # *strongly* suggested for sharing
    version='0.1',
    description='An example of a python package from pre-existing code',
    include_package_data=True,
    package_data={'cosmic_kite':['data/*']}
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)
