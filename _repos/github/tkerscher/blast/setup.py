from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='astroblast',
    version='1.0.2',
    description='Synchtrotron peak estimator for blazars',
    long_description='Estimates the synchrotron peak of blazars with prediction interval given a sed as produced by the VOUBlazar tool. Uses an ensemble of neural networks powered by pytorch.',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy'
    ],
    keywords='blazar synchrotron',
    url='https://github.com/tkerscher/blast',
    author='Tobias Kerscher',
    author_email='88444139+tkerscher@users.noreply.github.com',
    license='MIT',
    packages=['blast'],
    install_requires=[
        'numpy',
        'scipy',
        'torch',
        'tqdm'
    ],
    entry_points={
        'console_scripts': ['blast=blast.blast:main']
    },
    include_package_data=True,
    zip_safe=False)
