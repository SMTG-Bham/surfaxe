
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

Example: :mod:`surfaxe-parsefols --hkl 0,0,1 -b 8.83099`. 

