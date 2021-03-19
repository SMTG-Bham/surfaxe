[![Build Status](https://travis-ci.com/SMTG-UCL/surfaxe.svg?branch=master)](https://travis-ci.com/SMTG-UCL/surfaxe)  [![Documentation Status](https://readthedocs.org/projects/surfaxe/badge/?version=latest)](https://surfaxe.readthedocs.io/en/latest/?badge=stable) 
 
<img src='example_data/figures/surfaxe_header_v2.png' alt='surfaxe logo header' width='600'/>

Calculating the surface properties of crystals from first principles typically introduces several extra parameters including slab thickness, vacuum size, Miller index, surface termination and more. 
These factors all influence key properties of interest, making it a challenge to carry out simulations repeatably and draw reliable conclusions.
Surfaxe is a Python package for automating and simplifying density functional theory (DFT) calculations of surface properties, as well as providing analytical tools for bulk and surface calculations. 

The code is organised according to the best-practice workflow below:

<img src="example_data/figures/surfaxe_workflow_v4.png" alt="drawing" width="600"/>

The main features include:

1. **Slab generation:** Automated generation of surface slabs from command line.

  * All unique zero-dipole symmetric terminations of slabs are cleaved from a bulk structure.
  * Slabs can be organised into separate folders, optionally with all the required input files to run each calculation.

2. **Raw data processing:** Extracting data from convergence tests.

  * Parsing the convergence testing folders created using the slab generation scripts.
  * Producing dataframes and csv files of summary data.
  * Plotting scripts visualising convergence with respect to slab and vacuum thickness.

3. **Analysis:** Various scripts for surface and bulk calculations.

  * Calculation of planar and macroscopic average of the electrostatic potential through the slab to determine absolute electron energies (ionisation potential, electron affinity).
  * Nearest neighbour atom determination and bond distance analysis (useful for geometry relaxation convergence checks). 

Surfaxe primarily supports the [VASP](https://www.vasp.at/) DFT code, however most of the `generation` module is code-agnostic. In the future we would like to add support for more periodic codes in the other modules.

## Example outputs

**Analysis of average bond lengths as a function of slab thickness**
![Bond analysis example](example_data/figures/bond_analysis_plot.png)

**Surface energy convergence checks with respect to vacuum and slab thickness**
![Surface energy convergence example](example_data/figures/110_surface_energy.png)

See the [tutorials directory](tutorials/) for more examples.

## Installation

Surfaxe is a Python 3 package and requires pymatgen and other standard scientific Python packages.

Recommended installation is to git clone and install with `pip`:

```sh
git clone https://github.com/brlec/surfaxe.git
cd surfaxe
pip install --user .
```

 For development work, `--user` can be replaced with `-e`, which creates links to the source folder so any changes to the code are reflected on the path.

## Usage

### Quick start

Surfaxe can be used via the command line and via Python API. [The docs](https://surfaxe.readthedocs.io/en/latest/) include information on both, and the built-in `-h` option is available in the command line interface for each of the scripts.

We recommend starting off by looking at the dedicated [tutorials](https://github.com/SMTG-UCL/surfaxe/tree/master/tutorials). These Jupyter notebooks will guide you through most of the functionality of the package. 

The tutorials can also be run interactively on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SMTG-UCL/surfaxe/HEAD?filepath=tutorials)

### Command line interface
The scripts can be separated into four modules that follow a typical surfaces workflow; these are generation, convergence, analysis and data, with added plotting functionality. The vast majority of `surfaxe` functionality is available via the command line interface, but the Python API allows for more flexibility.

Generation:

* `surfaxe-getall`: Generates all unique symmetric zero-dipole surface slabs up to a maximum Miller index specified in any format supported by pymatgen. Optionally provides all DFT calculation input files.
* `surfaxe-gethkl`: Generates all unique symmetric zero-dipole surface slabs for one specified Miller index in any format supported by pymatgen. Optionally provides all DFT calculation input files.

Convergence:

* `surfaxe-parsefols`: Parses the convergence folders set up with `surfaxe-getall` or `surfaxe-gethkl` where calculations were run with VASP. Extracts the relevant data and plots convergence graphs of the variation of surface energy and energy per atom with respect to slab and vacuum thickness.

Analysis:

* `surfaxe-bonds`: Parses the structure, looking for bonds between specified atoms.
* `surfaxe-simplenn`: Predicts the coordination environment of atoms for simple structures.
* `surfaxe-complexnn`: Predcts the coordination environment of atoms in more complex structures where the default prediction algorithm fails.
* `surfaxe-potential`: Calculates and optionally plots the planar and macroscopic potential of the slab along c-axis.
* `surfaxe-cartdisp`: Calculates the Cartesian displacements of atoms during relaxation from intial and final structures.

Plotting: 

* `surfaxe-plot-surfen` and `surfaxe-plot-enatom`: Plot the surface energy and energy per atom based on data from `surfaxe-parsefols` with individual customisability
* `surfaxe-plot-bonds`: Plots the bond distance with respect to fractional coordinate, based on `surfaxe-bonds`
* `surfaxe-plot-potential`: Plots the planar and macroscopic potential based on data already analysed with `surfaxe-potential`

Data:

* `surfaxe-core`: Collects the core energy level from the middle of a surface slab, based on supplied bulk core atom and the list of its nearest neighbours.
* `surfaxe-vacuum`: Collects the vacuum potential level from a VASP LOCPOT file.

Pymatgen issues warnings whenever the hash in a VASP POTCAR present does not match the one in their library. Only one warning of the same type will be issued by default. All warnings can be suppressed completely by adding the following to your script:

```python
import warnings
warnings.filterwarnings('ignore')
```

## Development notes

### Bugs, features and questions

Please use the Issue Tracker to report bugs or request features in the first instance.

```sh
git clone https://github.com/SMTG-UCL/surfaxe.git
cd surfaxe
pip install --user .
```

Contributions to interface with this package are most welcome. Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

* Code style should comply with [PEP8](http://www.python.org/dev/peps/pep-0008) where possible. [Google's house style](https://google.github.io/styleguide/pyguide.html)
is also helpful, including a good model for docstrings.
* Please use comments liberally!
* Add tests wherever possible, and use the test suite to check if everything still works.

### Tests

Unit tests are in the `tests` directory and can be run from the top directory using `pytest`. Please run these tests whenever submitting bug fix pull requests and include new tests for new features as appropriate.

We also use CI testing on [travis](https://travis-ci.com/github/SMTG-UCL/surfaxe). 

## License and how to cite

Surfaxe is free to use under the MIT License. Please cite [add link] if you use it in your research.

## Dependencies

Surfaxe relies primarily on [Pymatgen](pymatgen.org) for manipulating crystal structures and interfacing with the DFT codes.

## Detailed requirements

Surfaxe is compatible with Python 3.6+ and requires the following packages:

* [Pymatgen](https://pymatgen.org/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [Numpy](https://numpy.org/)

## Contributors

List of contributors:

- Katarina Brlec (@brlec)
- Daniel Davies (@dandavies99)
- David Scanlon (@scanlond)

Acknowledgements:

- Surfaxe has benefited from useful discussions with Seán Kavanagh (@kavanase), Graeme Watson (@wantsong), Luisa Herring-Rodriguez (@zccalgh), Christopher Savory (@cnsavory), Bonan Zhu (@zhubonan) and Maud Einhorn (@maudeinhorn). 
- Thanks to Keith Butler (@keeeto) for providing a starting point and validation examples for calculating the planar average electrostatic potential via the [macrodensity](https://github.com/WMD-group/MacroDensity) package.  
