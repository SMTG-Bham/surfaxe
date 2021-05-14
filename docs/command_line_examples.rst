
Command Line Interface (CLI)
============================

While :mod:`surfaxe` has a full python API (see our tutorials page for example usage) it also has an
intuitive command line interface (CLI). Below are some simple examples of what you can do with
:mod:`surfaxe` at the command line. 

You can get a full list of accepted flags and what they do for each command using 
the :mod:`--help` or :mod:`-h` flag, e.g.:

.. code:: 

    $ surfaxe-bonds -h

    > usage: surfaxe-bonds [-h] [-s STRUCTURE] [-b BOND [BOND ...]]
                     [--oxi-list OX_STATES_LIST [OX_STATES_LIST ...]]
                     [--oxi-dict OX_STATES_DICT] [--no-csv]
                     [--csv-fname CSV_FNAME] [--no-plot]
                     [--plt-fname PLT_FNAME] [-c COLOR] [--width WIDTH]
                     [--height HEIGHT] [--dpi DPI] [--yaml]

    Parses the structure looking for bonds between atoms. Check the validity of
    the nearest neighbour method on the bulk structure before using it on slabs.

    optional arguments:
    -h, --help            show this help message and exit
    -s STRUCTURE, --structure STRUCTURE
                          Filename of structure file in any format supported by
                          pymatgen (default: POSCAR
    -b BOND [BOND ...], --bond BOND [BOND ...]
                          List of elements e.g. Ti O for a Ti-O bond
    --oxi-list OX_STATES_LIST [OX_STATES_LIST ...]
                          Add oxidation states to the structure as a list 
                          e.g. 3 3 -2 -2 -2
    --oxi-dict OX_STATES_DICT
                          Add oxidation states to the structure as a dictionary
                          e.g. Fe:3,O:-2
    --no-csv              Prints data to terminal
    --csv-fname CSV_FNAME
                          Filename of the csv file (default: bond_analysis.csv)
    --no-plot             Turns off plotting 
    --plt-fname PLT_FNAME
                          Filename of the plot (default: bond_analysis.png)
    -c COLOR, --color COLOR
                          Color of the marker in any format supported by mpl
                          e.g. "#eeefff" hex colours starting with # need to be
                          surrounded with quotation marks
    --width WIDTH         Width of the figure in inches (default: 6)
    --height HEIGHT       Height of the figure in inches (default: 5)
    --dpi DPI             Dots per inch (default: 300)
    --yaml YAML           Read all args from a yaml config file. Completely 
                          overrides any other flags set

The behaviour of default parameters of the functions is extensively documented in 
the :mod:`surfaxe` python package section of the docs. 

=======================
Pre-processing commands
=======================

**surfaxe-generate-slabs**: Generates all unique slabs with specific Miller indices or 
up to a maximum Miller index for a set of slab and vacuum thicknesses. 

Example: :mod:`surfaxe-generate -s bulk_structure.cif --hkl 1,1,0 -t 20 40 -v 20 40 -f` generates
all slabs for the (1,1,0) direction for minimum slab and vacuum thicknesses of 20 Å and 40 Å. 
The :mod:`-f` option organises these into subdirectories with all required VASP input 
files required to run singleshot calculations uisng default settings. It includes all combinations 
for zero-dipole terminations with inversion symmetry. 
The directory structure produced is:

.. code::

    100/              <-- Miller index
      ├── 20_20_0/    <-- slab-thickness_vacuum-thickness_termination-number
      ├── 20_40_0/   
      ├── 40_20_0/
      └── 40_40_0/
        ├── POSCAR    <-- VASP files 
        ├── INCAR
        ├── POTCAR
        └── KPOINTS

*Note: The hkl flag must be comma-separated with no spaces and the list of thicknesses and 
vacuums must be space-separated.*

*Note: To use the :mod:`-f` option you must first set up the 
`pymatgen POTCAR environment <https://pymatgen.org/installation.html#potcar-setup>`_.* 

Similarly, to above the script can be modified to consider multiple Miller indices. 

Example: :mod:`surfaxe-generate -s bulk_structure.cif --hkl 1,1,0 1,1,1 -t 20 40 -v 20 40 -f` 
generates all (1,1,0) and (1,1,1) slabs with minimum slab and vacuum thicknesses of 20 Å and 40 Å. 

*Note: h,k,l are comma-separated with no spaces, while the two (or more) Miller indices are space-separated.*

Lastly, a maximum hkl value can be supplied as an integer so that the script finds all 
zero-dipole slabs up to that maximum Miller index. 

Example: :mod:`surfaxe-generate -s SnO2.cif --hkl 2 -t 20 40 -v 30` generates all slabs with Miller 
indices up to a maximum value of 2, with minimum slab thicknesses of 20 Å and of 40 Å, and 
minimum vacuum thickness of 30 Å. 

========================
Post-processing commands
========================

**surfaxe-parsefols**: Parses data produced by electronic structure codes once calculations
have been run in then directory structures produced by one of the pre-processing commands. 

Example: :mod:`surfaxe-parsefols --hkl 0,0,1 -b 8.83099` saves a csv file of surface energies
and energies per atom for each slab-vacuum combination, as well as plots for each. See the 
Tutorials directory for examples. 

**surfaxe-plot-surfen** and **surfaxe-plot-enatom** can be used to customise the surface 
energy and energy per atom plots independetnly based on the data already collated 
with **surfaxe-parsefols**. 

=================
Analysis commands
=================

**surfaxe-potential**: Reads the local electrostatic potential file and plots the planar 
and macroscopic averages normal to the surface. Currently only the VASP LOCPOT 
file is supported as input. 

Example: :mod:`surfaxe-potential -l LOCPOT -v 11.5` produces a plot assuming a lattice vector of 
11.5 Angstroms and saves the plot data to a csv file. 

**surfaxe-bonds**: Analyse bonding in the structure using Pymatgen's local_env module.
Average bond lengths for each pair of species of interest can be plotted as a function 
of c lattice vector (normal to the slab surface). This can be useful for checking whether
the center of the slab has converged, where bond distances should be bulk-like. 

Example: :mod:`surfaxe-bonds -s CONTCAR -b Sn O` plots the average Sn-O bond length from the 
VASP output structure file. A csv file of the data plotted is also produced. 

**surfaxe-plot-potential** and **surfaxe-plot-bonds** can be used to generate the  
plots based on the data collated with **surfaxe-potential** and **surfaxe-bonds**, 
allowing customisation of plots without having to re-analyse the data. All plotting 
functionality is accessible through the main functions as well. 

**surfaxe-simplenn** and **surfaxe-complexnn**: Analyse the bonding in the slab, again using Pymatgen 
functions. *simplenn* is faster, but less reliable for systems with more complex bonding.
*complexnn* is more robust but requires a dictionary of cutoff bond lengths to be supplied
for each pair of species. See the analysis tutorial for further explanation. 

Example: :mod:`surfaxe-complexnn -s CONTCAR_bivo4 -b Bi3+ O2- 2.46 V5+ O2- 1.73` will 
analyse the coordination of atoms in this BiVO4 slab and save them to a csv file. 

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

================
YAML input files
================

Most CLI commands allow use of YAML input files containing all the arguments which cannot be 
used in conjunction with other command line argument flags. This is done by specifying 
the :mod:`--yaml` flag which overrides any other flags set in command line by loading the 
:mod:`surfaxe_config.yaml` file.

Sample YAML input files for each of the functions, with defaults and comments are in 
the :mod:`surfaxe/cli/templates` folder. 
All :mod:`**kwargs` of the main function can be passed in the YAML file.  

Example: Generation of (1,0,1) CdTe slabs could easily customised so that all VASP 
input files are created with specific INCAR tags using the following config.yaml file: 

.. code-block:: yaml 

    structure: CdTe.cif 
    hkl: (1,0,1) 
    thicknesses: [20, 40] 
    vacuums: [20, 40] 
    make_fols: True 
    make_files: True 
    max_size: 500 
    center_slab: True 
    ox_states: 
      Cd: 2
      Te: -2
    fmt: poscar 
    name: POSCAR 
    config_dict: PBE_config.json 
    user_incar_settings: 
      ENCUT: 460
      KPAR: 3
    user_kpoints_settings: 
      reciprocal_density: 35

The slabs would then be generated using :mod:`surfaxe-gethkl --yaml config.yaml`