#!/usr/bin/env python

__author__ = "Katarina Brlec"
__copyright__ = "Copyright Katarina Brlec, Daniel W. Davies (2020)"
__version__ = "0.1"
__maintainer__ = "Katarina Brlec"
__email__ = "katarina.brlec.15@ucl.ac.uk"
__date__ = "July 12 2020"

from setuptools import setup, Extension
import os

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
      name='surfaxe',
      version='1.0',
      description='Dealing with slabs for first principles calculations of surfaces',
      long_description=open(os.path.join(module_dir, 'README.md')).read(),
      long_description_content_type='text/markdown',
      url='https://github.com/SMTG-group/',
      author='Katarina Brlec',
      author_email='katarina.brlec.15@ucl.ac.uk',
      license='MIT',
      packages=['surfaxe'],
      zip_safe=False,
      install_requires=['scipy', 'numpy', 'spglib', 'pymatgen','pandas','scikit-learn'],
      python_requires='>=3.7',
      classifiers=[
        'Programming Language :: Python',
        'Development Status :: Development',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering'
      ],
      entry_points={
        'console_scripts': [
          'surfaxe-parse-energy = surfaxe.cli.parseenergy:main', 
          'surfaxe-parse-structures = surfaxe.cli.parsestructures:main', 
          'surfaxe-bonds = surfaxe.cli.bonds:main', 
          'surfaxe-simplenn = surfaxe.cli.simplenn:main', 
          'surfaxe-complexnn = surfaxe.cli.complexnn:main', 
          'surfaxe-potential = surfaxe.cli.potential:main', 
          'surfaxe-cartdisp = surfaxe.cli.cartdisp:main', 
          'surfaxe-core = surfaxe.cli.core:main', 
          'surfaxe-vacuum = surfaxe.cli.vacuum:main', 
          'surfaxe-plot-surfen = surfaxe.cli.plotsurfen:main', 
          'surfaxe-plot-enatom = surfaxe.cli.plotenatom:main', 
          'surfaxe-plot-bonds = surfaxe.cli.plotbonds:main', 
          'surfaxe-plot-potential = surfaxe.cli.plotpotential:main', 
          'surfaxe-generate = surfaxe.cli.generateslabs:main'
        ]
      }
    )
