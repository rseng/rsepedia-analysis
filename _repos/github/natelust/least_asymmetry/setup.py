# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os

cpp_args = ['-std=c++11', '-stdlib=libc++', '-mmacosx-version-min=10.7']

ext_modules = [
    Extension(
        'least_asymmetry.make_asym',
        ['least_asymmetry/make_asym.cc'],
        include_dirs=['least_asymmetry/include'],
        language='c++',
    ),
]


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.
    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }
    if sys.platform == 'darwin':
        c_opts['unix'].append('-stdlib=libc++')
        c_opts['unix'].append('-mmacosx-version-min=10.7')

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
        if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('least_asymmetry/include/')

setup(
    name='least_asymmetry',
    version='0.40',
    author='Nate Lust',
    author_email='nlust@astro.princeton.edu',
    description='A module for calculating centers though least asymmetry',
    packages=['least_asymmetry'],
    package_data={"": extra_files,},
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
    url='https://github.com/natelust/least_asymmetry',
    download_url='https://github.com/natelust/least_asymmetry/archive/0.40.tar.gz',
    keywords=['astronomy', 'centering', 'astrophysics', 'centroid', 'least asymmetry'],
    classifiers=[],
)
