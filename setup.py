# coding: utf-8
# Distributed under the terms of the MIT License.

import os
from glob import glob
from typing import Dict

from setuptools import find_packages, setup

with open('koopmans/__init__.py', 'r') as f:
    [versionline] = [line for line in f.readlines() if 'version' in line]
    version = versionline.split('=')[-1].strip().strip("'")

with open('requirements/requirements.txt', 'r') as f:
    requirements = [line.strip() for line in f.readlines()]

with open('requirements/pip_requirements.txt', 'r') as f:
    requirements += [line.strip() for line in f.readlines()]

with open('requirements/test_requirements.txt', 'r') as f:
    requirements += [line.strip() for line in f.readlines()]

extra_requirements: Dict[str, str] = dict(all=[])
req_files = glob('requirements/*.txt')
for _file in req_files:
    if _file not in ['requirements/requirements.txt', 'requirements/pip_requirements.txt']:
        with open(_file, 'r') as f:
            subreq = _file.split('/')[-1].split('_')[0]
            extra_requirements[subreq] = [line.strip() for line in f.readlines()]
            extra_requirements['all'] += extra_requirements[subreq]

with open("README.rst", "r") as f:
    long_description = f.read()

setup(name='koopmans',
      version=version,
      description='Koopmans spectral functional calculations with python and Quantum ESPRESSO',
      long_description=long_description,
      url='https://github.com/epfl-theos/koopmans',
      author='Edward Linscott',
      author_email='edwardlinscott@gmail.com',
      maintainer='Edward Linscott',
      maintainer_email='edwardlinscott@gmail.com',
      license='GPL',
      packages=find_packages() + find_packages('ase'),
      package_dir={'koopmans': 'koopmans', 'ase': 'ase/ase'},
      python_requires='>=3.7',
      install_requires=requirements,
      scripts=[s for s in glob('bin/*') if s[-2:] != '.x'],
      test_suite='tests',
      include_package_data=True,
      setup_requires=["setuptools>=42"],
      extras_require=extra_requirements,
      classifiers=[
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.7",
          "Programming Language :: Python :: 3.8",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Topic :: Scientific/Engineering",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Scientific/Engineering :: Physics"
      ],
      entry_points={'console_scripts': ['koopmans=koopmans.cli.main:main']},
      zip_safe=False)
