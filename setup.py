#!/usr/bin/env python
from setuptools import setup

setup(name='nonpher',
      version='1.0.1',
      description='Nonpher: computational method for design of hard-to-synthesize structures',
      url='https://github.com/lich-uct/nonpher',
      author='Milan Vorsilak',
      license='GPL-3.0',
      packages=['nonpher',], # 'molpher-lib', 'rdkit'],
      zip_safe=False,
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'nonpher = nonpher.nonpher:main',
          ],
      },
     )
