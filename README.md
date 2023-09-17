[![Build status](https://github.com/smtg-bham/surfaxe/actions/workflows/tests.yml/badge.svg)](https://github.com/SMTG-bham/surfaxe/actions)
[![Documentation Status](https://readthedocs.org/projects/surfaxe/badge/?version=latest)](https://surfaxe.readthedocs.io/en/latest/?badge=stable) 
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03171/status.svg)](https://doi.org/10.21105/joss.03171)

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
![Surface energy convergence example](example_data/figures/010_surface_energy.png)

See the [tutorials directory](tutorials/) for more examples.

## Installation

Surfaxe is a Python 3 package and requires pymatgen and other standard scientific Python packages.

Recommended installation is to git clone and install with `pip` in a stable virtual environment:

```sh
git clone https://github.com/SMTG-bham/surfaxe.git
cd surfaxe
pip install -e .
```

The `-e` option creates links to the source folder so any changes to the code are reflected on the path.

For the code to generate VASP input files along with the surface slabs, POTCARs need to be [set up with pymatgen](https://pymatgen.org/installation.html#potcar-setup).

## Usage

### Quick start

Surfaxe can be used via the command line and via Python API. [The docs](https://surfaxe.readthedocs.io/en/latest/) include information on both, and the built-in `-h` option is available in the command line interface for each of the scripts.

We recommend starting off by looking at the dedicated [tutorials](https://github.com/SMTG-bham/surfaxe/tree/master/tutorials). These Jupyter notebooks will guide you through most of the functionality of the package.

The tutorials can also be run interactively on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SMTG-bham/surfaxe/HEAD?filepath=tutorials)

### Command line interface

The scripts can be separated into four modules that follow a typical surfaces workflow; these are generation, convergence, analysis and data, with added plotting functionality. The vast majority of `surfaxe` functionality is available via the command line interface, but the python API allows for more flexibility.

Generation:

* `surfaxe-generate`: Generates all unique symmetric zero-dipole surface slabs for one or more specified Miller indices or up to a maximum Miller index specified, in any format supported by pymatgen. Optionally provides all VASP input files.

Convergence:

* `surfaxe-parse-energies`: Extracts the relevant data from the convergence folders set up with `surfaxe-generate` where calculations were run with VASP. Plots convergence graphs of the variation of surface energy with respect to slab and vacuum thickness. Can optionally parse core atom and vacuum energies as well.
* `surfaxe-parse-structures`: Collects the structures' metadata into a json file and optionally performs bond analysis as in `surfaxe-bonds`.

Analysis:

* `surfaxe-bonds`: Parses the structure, looking for bonds between specified atoms.
* `surfaxe-simplenn`: Predicts the coordination environment of atoms for simple structures.
* `surfaxe-complexnn`: Predcts the coordination environment of atoms in more complex structures where the default prediction algorithm fails.
* `surfaxe-potential`: Calculates the planar potential of the slab along c-axis, the gradient of the planar potential and optionally macroscopic potential.
* `surfaxe-surface-dipole`: Returns the surface dipole needed for macroscopic ionisation potential calculation
* `surfaxe-cartdisp`: Calculates the Cartesian displacements of atoms during relaxation from intial and final structures.

Plotting:

* `surfaxe-plot-surfen` and `surfaxe-plot-enatom`: Plot the surface energy and energy per atom based on data from `surfaxe-parsefols` with individual customisability
* `surfaxe-plot-bonds`: Plots the bond distance with respect to fractional coordinate, based on `surfaxe-bonds`
* `surfaxe-plot-potential`: Plots the planar (and macroscopic) potential based on data already analysed with `surfaxe-potential`

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

Contributions to interface with this package are most welcome. Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

* Code style should comply with [PEP8](http://www.python.org/dev/peps/pep-0008) where possible. [Google's house style](https://google.github.io/styleguide/pyguide.html)
is also helpful, including a good model for docstrings.
* Please use comments liberally!
* Add tests wherever possible, and use the test suite to check if everything still works.

### Tests

Unit tests are in the `tests` directory and can be run from the top directory using `pytest`. Please run these tests whenever submitting bug fix pull requests and include new tests for new features as appropriate.

We also use CI build and testing using [GitHub Actions](https://github.com/SMTG-bham/surfaxe/actions).

## License and how to cite

Surfaxe is free to use under the MIT License. If you use it in your research, please cite
> K. Brlec, D. W. Davies and D. O. Scanlon, *Surfaxe: Systematic surface calculations.* Journal of Open Source Software, 6(61), 3171, (2021) [DOI: 10.21105/joss.03171](https://joss.theoj.org/papers/10.21105/joss.03171)

## Dependencies

Surfaxe relies primarily on [Pymatgen](pymatgen.org) for manipulating crystal structures and interfacing with the DFT codes.

## Detailed requirements

Surfaxe is compatible with python 3.7+ and requires the following packages:

* [Pymatgen](https://pymatgen.org/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [Numpy](https://numpy.org/)
* [Scikit-learn](https://scikit-learn.org/stable/)

## Contributors

List of contributors:

* Katarina Brlec (@brlec)
* Daniel Davies (@dandavies99)
* David Scanlon (@scanlond)

Acknowledgements:

* Surfaxe has benefited from useful discussions with Adam Jackson (@ajjackson), Seán Kavanagh (@kavanase), Graeme Watson (@wantsong), Luisa Herring-Rodriguez (@zccalgh), Christopher Savory (@cnsavory), Bonan Zhu (@zhubonan) and Maud Einhorn (@maudeinhorn).
* Thanks to Keith Butler (@keeeto) for providing a starting point and validation examples for calculating the planar average electrostatic potential via the [macrodensity](https://github.com/WMD-group/MacroDensity) package.
* We thank @eihernan, @pzarabadip and @danielskatz for taking the time to review the code and offering valuable suggestions for improvements.
