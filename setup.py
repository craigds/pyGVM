#!/usr/bin/env python
from distutils.core import setup

version_tuple = __import__('pygvm').VERSION
version = ".".join([str(v) for v in version_tuple])

setup(
    name='pygvm',
    description='''Simplified port of GVM for python. See http://www.tomgibara.com/clustering/fast-spatial/java-library for the original.''',
    version=version,
    author='Craig de Stigter',
    author_email='craig.ds@gmail.com',
    url='http://github.com/craigds/pygvm',
    packages=['pygvm'],
    classifiers = [
        "Development Status :: 4 - Beta",
        "Environment :: Web Environment",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Utilities",
    ]
)
