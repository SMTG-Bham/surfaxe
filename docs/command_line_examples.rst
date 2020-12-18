
Command Line Interface (CLI)
============================

While :mod:`surfaxe` has a full python API (see our tutorials for example usage) it also has an
intuitive command line interface (CLI). Below are some simple examples of what you can do with
:mod:`surfaxe` at the command line. 

You can get a full list of accepted flags and what they do for each command using 
the :mod:`--help` or :mod:`-h` flag, e.g.:

.. code:: 

    $ surfaxe-bonds -h

    > usage: surfaxe-bonds [-h] [-s STRUCTURE] [-b LIST_OF_BONDS]
                     [--oxstates-list OX_STATES_LIST]
                     [--oxstates-dict OX_STATES_DICT] [--no-csv]
                     [--csv-fname CSV_FNAME] [--no-plot]
                     [--plt-fname PLT_FNAME] [--dpi DPI] [--yaml]

    Parses the structure looking for bonds between atoms. Check the validity of
    the nearest neighbour method on the bulk structure before using it on slabs.

    optional arguments:
    -h, --help            show this help message and exit
    -s STRUCTURE, --structure STRUCTURE
                            Filename of structure file in any format supported by
                            pymatgen
    -b LIST_OF_BONDS, --bonds LIST_OF_BONDS
                            List of bonds as lists to compare in any order (e.g.
                            Y-O,Ti-S
    --oxstates-list OX_STATES_LIST
                            Add oxidation states to the structure as a list.
    --oxstates-dict OX_STATES_DICT
                            Add oxidation states to the structure as a dictionary
                            e.g. "Fe:3,O:-2"
    --no-csv              Turns off saving data to csv file
    --csv-fname CSV_FNAME
                            Filename of the csv file (default: bond_analysis.csv)
    --no-plot             Turns off plotting the bond lengths
    --plt-fname PLT_FNAME
                            Filename of the plot
    --dpi DPI             Dots per inch
    --yaml                Read optional args from surfaxe_config.yaml file.

=======================
Pre-processing commands
=======================

**surfaxe-gethkl**: Generates all unique slabs with a specific Miller index for a set of 
slab and vacuum thicknesses. 

Example: :mod:`surfaxe-gethkl -s bulk_structure.cif --hkl 1,0,0 -t 20 40 -v 20 40 -f` generates
all slabs for the (1,0,0) direction for minimum slab and vacuum thicnesses of 20 and 40. 
The :mod:`-f` option organises these into subdirectories with all required VASP input 
files required to run singleshot calculations uisng default settings. It includes all combinations 
for zero-dipole, symmetric terminations.
The directory structure produced is:

.. code::

    100/              <-- miller index
      ├── 20_20_0/    <-- slab-thickness_vacuum-thickness_termination-number
      ├── 20_40_0/   
      ├── 40_20_0/
      └── 40_40_0/
        ├── POSCAR    <-- vasp files 
        ├── INCAR
        ├── POTCAR
        └── KPOINTS

*Note: The hkl flag must be comma-separated and the list of thicknesses and vacuums must 
be space-separated.*


**surfaxe-getall**: Similar to above but considers multiple Miller indices. A maximum hkl value must be 
supplied as an integer.

Example: :mod:`surfaxe-getall -s SnO2.cif --hkl 1 -t 20 40 -v 30` generates all slabs with Miller indices 
up to a maximum value of 1, with minimum slab thicknesses of 20 and of 40, and minimum vacuum 
thickness of 30. 

========================
Post-processing commands
========================

**surfaxe-parsefols**: Parses data produced by electronic structure codes once calculations
have been run in then directory structures produced by one of the pre-processing commands. 

Example: :mod:`surfaxe-parsefols --hkl 0,0,1 -b 8.83099` saves a csv file of surface energies
and energies per atom for each slab-vacuum combination, as well as plots for each. See the 
Tutorials directory for examples. 

=================
Analysis commands
=================

**surfaxe-potential**: Reads the local electrostatic potential file and plots the planar and macroscopic
averages normal to the surface (inspired by PlanarAverage.py in  
`Keith Butler's Macrodensity code <https://www.github.com/WMD-group/macrodensity>`_. Currently
only the VASP LOCPOT file is supported as iput. 

Example: :mod:`surfaxe-potential -l LOCPOT -v 11.5` produces a plot assuming a lattice vector of 
11.5 Angstroms and saves the plot data to a csv file. 

**surfaxe-bonds**: Analyse bonding in the structure using Pymatgen's local_env module.
Average bond lengths for each pair of species of interest can be plotted as a function 
of c lattice vector (normal to the slab surface). This can be useful for checking whether
the center of the slab has converged, bulk-like bond distances. 

Example: :mod:`surfaxe-bonds -s CONTCAR -b Sn-O` plots the average Sn-O bond length from the 
VASP output structure file. A csv file of the data plotted is also produced. 

**surfaxe-simplenn** and **surfaxe-complexnn**: Analyse the bonding in the slab, again using Pymatgen 
functions. *simplenn* is faster, but less reliable for systems with more complex bonding.
*complexnn* is more robust but requires a dictionary of cutoff bond lengths to be supplied
for each pair of species. See the analysis tutorial for further explanation. 

=============
Data commands
=============

There are some simple convenience commands that can also be used to extract key values from
raw data files produced by solid state codes. Currently only commands relating to VASP output
files are included, which rely on the surfaxe :mod:`vasp_data` module. We hope to expand this
in the future. 

**surfaxe-vacuum** and **surfaxe-core** can be used to extract vacuum and core energies, respectively, 
that are needed to calculate absolute electron energies (ionisation potential and electron affinity). 
See the `Macrodensity <https://www.github.com/WMD-group/macrodensity>`_ tutorials for more information
on the steps needed to do this. 