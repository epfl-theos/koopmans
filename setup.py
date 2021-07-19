# coding: utf-8
# Distributed under the terms of the MIT License.

import os
from glob import glob
from setuptools import setup, find_packages
from koopmans import __version__

with open('requirements/requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f.readlines()]

with open('requirements/pip_requirements.txt', 'r') as f:
    requirements += [line.strip() for line in f.readlines()]

with open('requirements/test_requirements.txt', 'r') as f:
    requirements += [line.strip() for line in f.readlines()]

extra_requirements = dict(all=[])
req_files = glob('requirements/*.txt')
for _file in req_files:
    if _file not in ['requirements/requirements.txt', 'requirements/pip_requirements.txt']:
        with open(_file, 'r') as f:
            subreq = _file.split('/')[-1].split('_')[0]
            extra_requirements[subreq] = [line.strip() for line in f.readlines()]
            extra_requirements['all'] += extra_requirements[subreq]

with open("README.rst", "r") as f:
    long_description = f.read()

setup(name='python_KI',
      version=__version__,
      description='Koopmans spectral functional calculations with python and quantum espresso',
      long_description=long_description,
      url='https://github.com/elinscott/python_KI',
      author='Edward Linscott',
      author_email='edwardlinscott@gmail.com',
      maintainer='Edward Linscott',
      maintainer_email='edwardlinscott@gmail.com',
      license='MIT',
      packages=find_packages() + find_packages(where='./ase_koopmans'),
      package_dir={'': '.', 'ase': './ase_koopmans/ase'},
      python_requires='>=3.6',
      install_requires=requirements,
      scripts=glob('scripts/*'),
      test_suite='tests',
      include_package_data=True,
      setup_requires=["setuptools>=42"],
      extras_require=extra_requirements,
      classifiers=[
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.6",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Scientific/Engineering :: Physics"
      ],
      entry_points={'console_scripts': ['koopmans=koopmans.cli.main:main']},
      zip_safe=False)
