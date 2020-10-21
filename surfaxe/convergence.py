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
mpl.rcParams.update({'font.size': 14})

# surfaxe
from surfaxe.io import plot_enatom, plot_surfen

def slab_from_file(structure, hkl):
    """
    Reads in structure from the file and returns slab object.

    Args:
         structure (str): Structure file in any format supported by pymatgen.
         hkl (tuple): Miller index of the slab in the input file.

    Returns:
         Slab object
    """
    slab_input = Structure.from_file(structure)
    return Slab(slab_input.lattice,
                slab_input.species_and_occu,
                slab_input.frac_coords,
                hkl,
                Structure.from_sites(slab_input, to_unit_cell=True),
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=slab_input.site_properties)

def parse_fols(hkl=None, bulk_per_atom=None, path_to_fols=None, plt_enatom=True, 
plt_surfen=True, save_csv=True):
    """
    Parses the convergence folders to get the surface energy, total energy,
    energy per atom and time taken for each slab and vacuum thickness
    combination.

    Args:
        hkl (`tuple`, required): Miller index of the slab.
        bulk_per_atom (`float`, required): Bulk energy per atom from a converged 
            bulk calculation.
        plt_enatom (`bool`, optional): Plots the energy per atom. Defaults to 
            ``True``.
        plt_surfen (`bool`, optional): Plots the surface energy. Defaults to 
            ``True``.
        save_csv (`bool`, optional): Saves the csv. Defaults to ``True``.

    Returns:
        DataFrame 
    """
     # Check all neccessary input parameters are present 
    if not any ([hkl, bulk_per_atom]): 
        raise ValueError('One or more of the required arguments (bulk_per_atom,'
                         ' hkl) were not supplied.')
    
    df_list = []
    hkl_string = ''.join(map(str, hkl))

    if path_to_fols is not None:
        cwd = path_to_fols
    else: 
        cwd = os.getcwd()

    for root, fols, files in os.walk(cwd):
        for fol in fols:
            if not any([fol==root, fol=='.ipynb_checkpoints']):
                path = os.path.join(root, fol)
                vsp_path = '{}/vasprun.xml'.format(path)
                otc_path = '{}/OUTCAR'.format(path)

                # instantiate structure, slab, vasprun and outcar objects
                vsp = Vasprun(vsp_path)
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
    if plt_enatom: 
        plot_enatom(df, hkl=hkl)
    
    if plt_surfen: 
        plot_surfen(df, hkl=hkl)

    # Save the csv or return the dataframe
    if save_csv: 
        df.to_csv('{}_data.csv'.format(hkl_string), header=True, index=False)
    
    else: 
        return df


