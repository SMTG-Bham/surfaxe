# Pymatgen
from pymatgen.core.surface import Slab
from pymatgen import Structure, Specie, Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.analysis.local_env import CrystalNN

# Misc
import os
import math
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('once')

# Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# surfaxe
from surfaxe.io import plot_enatom, plot_surfen, slab_from_file
from surfaxe.vasp_data import vacuum, core_energy

def parse_fols(hkl, bulk_per_atom, path_to_fols=None, parse_core_energy=False,
core_atom=None, bulk_nn=None, parse_vacuum=False, plt_enatom=True, 
plt_enatom_fname='energy_per_atom.png', plt_surfen=True, 
plt_surfen_fname='surface_energy.png', save_csv=True, **kwargs):
    """
    Parses the convergence folders to get the surface energy, total energy,
    energy per atom, band gap and time taken for each slab and vacuum thickness
    combination. It can optionally parse vacuum and core level energies. 
    The convergence folders must be the only ones in the folder to 
    which the ``path`` leads. Folders must contain vasprun.xml and 
    OUTCAR files. 

    Args:
        hkl (`tuple`): Miller index of the slab.
        bulk_per_atom (`float`): Bulk energy per atom from a converged 
            bulk calculation in eV per atom.
        path_to_fols (`str`, optional): Relative path to the convergence folders. 
            Defaults to cwd
        parse_core_energy (`bool`, optional): If True the scripts attempts to 
            parse core energies from a supplied OUTCAR. Defaults to ``False``. 
        core_atom (`str`, optional): The symbol of atom the core state energy 
            level should be parsed from. Defaults to ``None``. 
        bulk_nn (`list`, optional): The symbols of the nearest neighbours of the 
            `core_atom`. Defaults to ``None``. 
        parse_vacuum (`bool`, optional): if ``True`` the script attempts 
            to parse LOCPOT using analysis.electrostatic_potential to use the 
            maximum value of planar potential as the vacuum energy level. 
            Defaults to ``True``. 
        plt_enatom (`bool`, optional): Plots the energy per atom. Defaults to 
            ``True``.
        plt_enatom_fname (`str`, optional): The name of the energy per atom plot. 
            Defaults to ``energy_per_atom.png``.
        plt_surfen (`bool`, optional): Plots the surface energy. Defaults to 
            ``True``.
        plt_surfen_fname (`str`, optional): The name of the surface energy plot.
            Defaults to ``surface_energy.png``
        save_csv (`bool`, optional): Saves the csv. Defaults to ``True``.

    Returns:
        DataFrame 
    """
    
    # Update kwargs for core energy 
    get_core_energy_kwargs = {'orbital': '1s', 'ox_states': None, 
    'nn_method': CrystalNN()}
    get_core_energy_kwargs.update(
        (k, kwargs[k]) for k in get_core_energy_kwargs.keys() & kwargs.keys()
    )

    df_list, electrostatic_list, core_energy_list = ([] for i in range(3))
    hkl_string = ''.join(map(str, hkl))

    if path_to_fols is None:
        cwd = os.getcwd()
    else: 
        cwd = path_to_fols

    for root, fols, files in os.walk(cwd):
        for fol in fols:
            if not any([fol==root, fol=='.ipynb_checkpoints']):
                path = os.path.join(root, fol)
                vsp_path = '{}/vasprun.xml'.format(path)
                otc_path = '{}/OUTCAR'.format(path)

                # instantiate structure, slab, vasprun and outcar objects
                vsp = Vasprun(vsp_path, parse_potcar_file=False)
                otc = Outcar(otc_path)
                slab = slab_from_file(vsp.final_structure, hkl)
                vsp_dict = vsp.as_dict()

                # extract the data
                otc_times = otc.run_stats

                # name of fol has to be ./slabthickness_vacthickness_index
                slab_vac_index = fol.split('_')

                df_list.append(
                    {'hkl_string': hkl_string, 
                    'hkl_tuple': hkl, 
                    'slab_thickness': slab_vac_index[0],
                    'vac_thickness': slab_vac_index[1],
                    'slab_index': slab_vac_index[2],
                    'atoms': vsp_dict['nsites'], 
                    'area': slab.surface_area, 
                    'bandgap': vsp_dict['output']['bandgap'],
                    'slab_energy': vsp_dict['output']['final_energy'],
                    'slab_per_atom': vsp_dict['output']['final_energy_per_atom'],
                    'time_taken': otc_times['Elapsed time (sec)']})

                if parse_vacuum: 
                        electrostatic_list.append(
                            vacuum(path)
                        )
                                                
                if parse_core_energy: 

                    core_energy_list.append(
                        core_energy(core_atom, bulk_nn, outcar=otc_path, 
                        structure=slab, **get_core_energy_kwargs)
                        ) 

    df = pd.DataFrame(df_list)
    df['surface_energy'] = (
        (df['slab_energy'] - bulk_per_atom * df['atoms'])/(2*df['area']) * 16.02
        ) 

    if electrostatic_list: 
        df['vacuum_potential'] = electrostatic_list
    
    if core_energy_list: 
        df['core_energy'] = core_energy_list

    # Plot energy per atom and surface energy
    plt_kwargs = {'time_taken': True, 'cmap': 'Wistia', 'dpi': 300, 
    'heatmap': False}
    plt_kwargs.update((k, kwargs[k]) for k in plt_kwargs.keys() & kwargs.keys())

    if plt_enatom: 
        plot_enatom(df, plt_fname=plt_enatom_fname, **plt_kwargs)
    
    if plt_surfen: 
        plot_surfen(df, plt_fname=plt_surfen_fname, **plt_kwargs)

    # Save the csv or return the dataframe
    if save_csv: 
        df.to_csv('{}_data.csv'.format(hkl_string), header=True, index=False)
    
    else: 
        return df


