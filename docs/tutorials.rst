Tutorials
=========

We recommend starting off by looking at the `dedicated tutorials. <https://github.com/SMTG-UCL/surfaxe/tree/master/tutorials>`_ 
These Jupyter notebooks will guide you through most of the functionality of the package. 

The tutorials can also be run interactively `on Binder. <https://mybinder.org/v2/gh/SMTG-UCL/surfaxe/HEAD?filepath=tutorials>`_


================================
Using configuration dictionaries
================================

One of the most powerful parts of surfaxe is its ability to make all VASP input 
files needed for convergence testing. To do so surfaxe makes use of configuration
dictionaries (config dicts for short). These are python dictionaries that contain 
information used to set up INCAR, KPOINTS and POTCAR files. 

For example, if we were interested in setting up a single shot PBEsol calculation 
on SnO2 slabs, we could set up the config dict as follows:

.. code-block:: python

    config_dict = {
    "INCAR": {
        "ALGO": "Normal",
        "EDIFF": 1e-06,
        "EDIFFG": -0.01,
        "ENCUT": 500,
        "GGA": "PS",
        "ISMEAR": 0,
        "ISYM": 2,
        "IWAVPR": 1,
        "LASPH": true,
        "LORBIT": 11,
        "LREAL": "auto",
        "NELM": 200,
        "NSW": 0,
        "PREC": "Accurate",
        "SIGMA": 0.02
    },
    "KPOINTS": {
        "reciprocal_density": 55
    },
    "POTCAR": {
        "Sn": "Sn_d", 
        "O" : "O"
    }
    }

Alternatively, one of the ready-made :mod:`surfaxe` config dicts (:mod:`PBEsol.json`, 
:mod:`PBEsol_relax.json`, :mod:`PBE.json`, :mod:`PBE_relax.json` or :mod:`HSE06.json`) 
can be used and further modified using :mod:`user_incar_settings`, 
:mod:`user_kpoints_settings` and :mod:`user_potcar_settings`. The :mod:`relax` config dicts 
contain additional parameters necessary for geometric relaxations of slabs. 
The POTCAR functional (i.e. PBE, PBE_54) can be chosen with :mod:`user_potcar_functional`.  

`Pymatgen documentation <https://pymatgen.org/pymatgen.io.vasp.sets.html#pymatgen.io.vasp.sets.DictSet>`_ 
covers exact behaviour of the :mod:`user_incar_settings`, :mod:`user_kpoints_settings` and :mod:`user_potcar_settings` 
and all additional keyword arguments that can be supplied to slab generation scripts.  
