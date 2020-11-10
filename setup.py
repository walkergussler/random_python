#!/usr/bin/env python
from setuptools import setup,find_packages

import versioneer

long_descrption='''\
'''

setup(
  name='phacelia',
  version=versioneer.get_version(),
  description='phacelia recent/chronic hcv caller',
  long_description=long_descrption,
  author='CDC OID/NCHHSTP/DVH Bioinformatics Division',
  packages=find_packages(),
  cmdclass=versioneer.get_cmdclass(),
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Programming Language :: Python :: 2.7'
  ],
  setup_requires=[
  ],
  install_requires=[
    'six',
    'pandas',
    'numpy',
    'scipy',
    'scikit-learn',
    'toolz'
  ],
  entry_points={
    'console_scripts': [
      'phacelia = phacelia.entry:main'
    ]
  },
  package_data={'phacelia' : [
    'data/PHACELIA.pkl',
    'data/dnaDAC_10.txt'
    ]
  },
  include_package_data=True,
  test_suite='tests',
)
