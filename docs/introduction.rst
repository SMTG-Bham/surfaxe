
Getting started 
===============

============
Introduction
============

The purpose of :mod:`surfaxe` is to simplify the calculation 
of surface properties of inorganic compounds using first-principles codes. 
The package includes pre-processing tools to cleave surfaces from the bulk, automatically 
set up calculation directories and facilitate convergence testing, as well as
post-processing tools to calcualte surface energies, electrostatic potentials, and produce
publication-qaulity plots. 

============
Installation
============

Installation will be via :mod:`pip` in the near future. For now, please clone the repository 
and install the latest stable version:

.. code:: bash

    git clone https://github.com/SMTG-UCL/surfaxe.git
    cd surfaxe
    pip install -e .

The :mod:`-e` is optional and will install the project in developer (editable) mode.

For the code to generate VASP input files along with the surface slabs, 
`pymatgen POTCAR environment <https://pymatgen.org/installation.html#potcar-setup>`
must be set up correctly. 

=====
Usage
=====

In general, there are two ways to use :mod:`surfaxe`: 
(i) at the command line or (ii) using python in scripts or notebooks. 

Take a look at our `command line interface (CLI) examples <command_line_examples.html>`_ for an overview
of the CLI tools available. 

There is full documentation for all modules (have a browse of the side bar on the left)  
if you would rather use the python API directly. 