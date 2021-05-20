---
title: 'Surfaxe: Systematic surface calculations'
tags:
  - Python
  - materials design
  - high-throughput screening
  - surfaces
authors:
  - name: Katarina Brlec
    orcid: 0000-0003-1485-1888
    affiliation: 1, 2
  - name: Daniel W. Davies
    orcid: 0000-0003-4094-5992
    affiliation: 1, 2 
  - name: David O. Scanlon
    orcid: 0000-0001-9174-8601
    affiliation: 1, 2, 3
affiliations:
 - name: Department of Chemistry, University College London, 20 Gordon Street, London WC1H 0AJ, United Kingdom
   index: 1
 - name: Thomas Young Centre, University College London, Gower Street, London WC1E 6BT, United Kingdom
   index: 2
 - name: Diamond Light Source Ltd., Diamond House, Harwell Science and Innovation Campus, Didcot, Oxfordshire OX11 0DE, UK
   index: 3
date: 01 January 2021
bibliography: paper.bib
---

# Summary 
Surface science is key to understanding the properties of a wide range of materials for energy applications, from catalysts to solar cells to battery components. Computational modelling based on quantum mechanics is often used to calculate surface properties of materials, which in turn determine their stability and performance. The maturity of these "first-principles" methods, coupled with the huge amount of computational power accessible today, means they can now be used predictively in high-throughput screening workflows to suggest new materials for specific applications before they are synthesised. The `surfaxe` package provides a framework for such screening workflows, automating each stage of the process. 

# Statement of need 
Accurate descriptions of electronic structure are needed to calculate surface properties including surface formation energies, adsorption energies and absolute electron energies (ionisation potential, electron affinity, and work function). However, surface calculations diverge significantly from the typical setup in periodic codes where the bulk crystal is described by repeat boundary conditions in all three dimensions. To reveal a surface, the bulk must be cleaved into a slab and periodicity reduced to just two dimensions. Additional parameters must then be considered such as variance in slab thickness, vacuum size, surface termination (dangling bonds), and net electrostatic dipole moments, all of which complicate the calculation workflow and make reliable determination of properties more difficult.

# Surfaxe
The aims of `surfaxe` are:

- To act as a framework for the automation of surface calculations, with particular emphasis on ensuring that properties are converged with respect to the additional parameters that are introduced compared to bulk calculations
- To increase the efficiency and reproducibility of surface calculations by automating the generation of input files and processing of output files for density functional theory (DFT) codes
- To provide a toolbox of intuitive analytical tools to calculate performance-critical materials properties and directly generate publication-quality plots

The code makes extensive use of existing Python Materials Genomics (pymatgen) [@pymatgen] surface modules with full functionality retained in `surfaxe`. As well as a fully flexible Python API, `surfaxe` has a lightweight command line interface.
The modularity of `surfaxe` closely follows a best-practice workflow for the calculation of surface properties, with key features including:

- Automatic cleaving of slabs from the bulk crystal and organising them into a directory structure with all necessary calculation input files -- generation module
- Analyses of atomic displacements and coordination environments, bond lengths and electrostatic potential through the slab (Figure 1a and b) -- analysis module
- Processing of raw DFT outputs to determine surface energy variation with slab and vacuum thickness (Figure 1c) -- convergence module
- Automatic extraction of surface energy, vacuum and core energy levels, along with the important calculation parameters -- data module

In addition to pymatgen, existing packages related to surface calculations include the Atomic Simulation Environment (ase) [@ase], which is a large materials informatics library, and smaller packages to aid with specific post-processing tasks: MacroDensity for plotting of potentials [@macrodensity], WullfPack for plotting of Wulff shapes [@wulffpack], and bapt for plotting band alignments [@bapt]. While these toolkits are extremely useful, `surfaxe` is distinct with its focus  on the rigorous convergence of properties, the enabling of reproducible workflows, and the production of processed datasets and plots at the command line. Lastly, `surfaxe` is built on the pymatgen ecosystem, so full integration with the workflow packages FireWorks [@fireworks] and AiiDA [@aiida] is possible for managing calculations on high-performance computing clusters. 

![Example analysis: a) average bond length, b) electrostatic potential as a function of lattice parameter perpendicular to the surface, and c) a typical surface energy convergence plot with respect to slab and vacuum thickness. \label{fig1}](figures/joss_fig1.png)

# Acknowledgements
The development of this code has benefited from useful discussions with Seán Kavanagh, Graeme W. Watson, Luisa Herring-Rodriguez, Christopher N. Savory, Bonan Zhu, and Maud Einhorn. KB, DWD, and DOS acknowledge support from the European Research Council, ERC, (Grant 758345).



# References
