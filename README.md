# Surfaxe (neé slabby-stabby)

Surfaxe is a python package for automating and simplifying surface calculations, as well as provide some analytical tools for bulk and surface calculations. I kinda gave up on writing this but at least it's a start of some sort.

The main features include:
1. Organised generation of surface slabs from command line.
  * All unique zero-dipole symmetric terminations of slabs are cleaved from the bulk using Pymatgen.
  * The surface slabs can be generated in separate labelled folders with or without all VASP input files.

2. Extracting data from surfaces convergence testing created with Slabaxe.
  * Parses the convergence testing folders, using Pymatgen to read VASP outputs.
  * Plotting scripts for surface energy and energy per atom convergence from the output files created with parsing module.

3. Various analysis scripts for surface and bulk calculations.
  * Electrostatic potential tool is based on MacroDensity [link]
  * Nearest neighbours and bond analysis scripts

Surfaxe primarily supports VASP, however we would like to add support for other solid-state codes in the future. Code contributions to interface with these packages are welcome. See the contributing section for information about reporting bugs or getting involved.

Surfaxe is free to use, however we ask you cite [add link] it if you use it in your research.

## Usage
Slabaxe can be used via the command line and python API. The manual, including tutorials and API documentation can be found here. Additionally, the built-in `- h` option is available in the command line interface for each of the scripts.

The Tutorial page contains a detailed guide to using each command.

For a preview of the functionality, see the Gallery.

Currently the provisional scripts are separated into three stages; creation, convergence and analysis.
Creation:
* `surfaxe-getall`: For generating all unique symmetric zero-dipole surface slabs up to a maximum Miller index specified.
* `surfaxe-getone`: For generating all unique symmetric zero-dipole surface slabs for one specified Miller index.

Convergence:
* `surfaxe-parse`: For parsing the convergence folders set up with `surfaxe-getall` or `surfaxe-getone`.
* `surfaxe-surfenplot`: For plotting the surface energy and time taken for calculation to complete for each unique termination of surface slabs in convergence tests.
* `surfaxe-enatomplot`: For plotting the energy per atom and time taken for calculation to complete for each unique termination of surface slabs in convergence tests.

Analysis:
* `surfaxe-bond`
* `surfaxe-bondplot`
* `surfaxe-simplenn`
* `surfaxe-complexnn`
* `surfaxe-electrostatic`

## Installation
Surfaxe is a Python 3 package and requires pymatgen and other standard scientific python dependencies. While Surfaxe is a relatively small package, the dependencies are about
I think this should work?
```
pip install --user surfaxe
```


### Developer Installation
The source code for Surfaxe can be cloned from the git repo. For development work this is the preferred way of installation as it creates links to the source folder so any changes to the code are reflected on the path.

```sh
git clone https://github.com/brlec/slabby-stabby.git
cd Surfaxe
pip install --user -e .
```
#### Tests
no such thing as tests, we blindly trust a 23yo gal with four months experience!

## License

## Detailed requirements
Slabaxe is compatible with python 3.something+ and requires the following packages:
* Pymatgen
* Numpy
* Pandas
* Matplotlib
