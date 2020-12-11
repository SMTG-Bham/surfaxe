# Pymatgen
from pymatgen.core.surface import Slab
from pymatgen import Structure, Specie, Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.outputs import Vasprun, Outcar

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

def parse_fols(hkl, bulk_per_atom, path_to_fols=None, plt_enatom=True, 
plt_surfen=True, save_csv=True, **kwargs):
    """
    Parses the convergence folders to get the surface energy, total energy,
    energy per atom and time taken for each slab and vacuum thickness
    combination. The convergence folders must be the only ones in the folder to 
    which the ``path`` leads. Folders need to contain vasprun.xml and 
    OUTCAR files. 

    Args:
        hkl (`tuple`): Miller index of the slab.
        bulk_per_atom (`float`): Bulk energy per atom from a converged 
            bulk calculation in eV per atom.
        path_to_fols (`str`, optional): Relative path to the convergence folders. 
            Defaults to cwd
        plt_enatom (`bool`, optional): Plots the energy per atom. Defaults to 
            ``True``.
        plt_surfen (`bool`, optional): Plots the surface energy. Defaults to 
            ``True``.
        save_csv (`bool`, optional): Saves the csv. Defaults to ``True``.

    Returns:
        DataFrame 
    """
   
    df_list = []
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
                slab = slab_from_file(vsp_path, hkl)
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
                    'slab_energy': vsp_dict['output']['final_energy'],
                    'slab_per_atom': vsp_dict['output']['final_energy_per_atom'],
                    'time_taken': otc_times['Elapsed time (sec)']})

    df = pd.DataFrame(df_list)
    df['surface_energy'] = (
        (df['slab_energy'] - bulk_per_atom * df['atoms'])/(2*df['area']) * 16.02
        ) 

    # Plot energy per atom and surface energy
    plt_kwargs = {'time_taken': True, 'cmap': 'Wistia', 'fmt': 'png', 'dpi': 300, 
    'heatmap': False}
    plt_kwargs.update(**kwargs)

    if plt_enatom: 
        plot_enatom(df, hkl=hkl, **plt_kwargs)
    
    if plt_surfen: 
        plot_surfen(df, hkl=hkl, **plt_kwargs)

    # Save the csv or return the dataframe
    if save_csv: 
        df.to_csv('{}_data.csv'.format(hkl_string), header=True, index=False)
    
    else: 
        return df


