
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
                            Add oxidation states to the structure asa dictionary
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

**surfaxe-getall**: Generates all unique slabs with specified maximum Miller index, minimum slab
and vacuum thicknesses. It includes all combinations for multiple zero dipole
symmetric terminations for the same Miller index.

Example: :mod:`surfaxe-getall -s SnO2.cif --hkl 1 -t 20,40 -v 30` generates all slabs
up to a max miller index of 1, with minimum slab thicknesses of 20 and of 40, and minimum vacuum 
thickness of 30. 