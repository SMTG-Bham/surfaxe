# Pymatgen  
from pymatgen.core.surface import Slab
from pymatgen import Structure, Specie, Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.outputs import Locpot, Outcar, Vasprun
from pymatgen.analysis.local_env import BrunnerNN_real, BrunnerNN_reciprocal,\
BrunnerNN_relative, CovalentBondNN, Critic2NN, CrystalNN, CutOffDictNN, EconNN,\
JmolNN, MinimumDistanceNN, MinimumOKeeffeNN, MinimumVIRENN, NearNeighbors, VoronoiNN

# Misc
import os
import pandas as pd 
import numpy as np 
import warnings 
from surfaxe.generation import custom_formatwarning

# surfaxe 
from surfaxe.convergence import slab_from_file
from surfaxe.generation import oxidation_states

def data_collection(bulk_per_atom=None, folders=None, hkl_dict=None, 
                    parse_folders=True, 
                    return_df=False, parse_core_energy=True, 
                    parse_electrostatic=True, **kwargs): 
    """
    
    Args:
    bulk_per_atom (float): bulk energy per atom from a converged bulk calculation
    folders (list): list of folders where the calculations are ? 
    hkl_dict (dict): dictionary of tuples of Miller indices and the folders the 
    calculations relating to them are in, e.g. {(1,-1,2): '1-12'}; default=None
    parse (bool): if true the script parses the hkl root folders; default=True
    return_df (bool): if true will return the dataframe, default=False
    parse_core_energy (bool): if True will attempt to parse core energies from a
    supplied OUTCAR 
    parse_electrostatic (bool): if True will attempt to parse LOCPOT using 
    analysis.electrostatic_potential() 
    Returns: 

    """
    # Check if hkl_dict is correctly set up 
    if hkl_dict: 
        for key, value in hkl_dict.items(): 
            if isinstance(key, tuple): 
                raise TypeError('The keys supplied to hkl_dict are not tuples.')
            if isinstance(value, str): 
                raise TypeError('The values supplied to hkl_dict are not strings.')

    # Get the Miller indices as tuples and strings from folders in root dir
    if parse_folders:
        if not hkl_dict: hkl_dict = {}
        for root, fols, files in os.walk('.'):
            for fol in fols: 
                if fol is not root and len(fol)==3 and fol.isdigit():
                    hkl_dict[tuple(map(int, fol))] = fol

    # For each miller index, check if the folders specified are there and 
    # parse them for data

    df_list, electrostatic_list, core_energy_list = ([] for i in range(3))

    for folder in folders: 
        for hkl_tuple, hkl_string in hkl_dict.items(): 
            for root, fols, files in os.walk('./{}'.format(hkl_string)): 
                for fol in fols: 
                    if folder == fol: 
                        path = os.path.join(root,fol)
                        vsp_path = '{}/vasprun.xml'.format(path)

                        vsp = Vasprun(vsp_path, **kwargs)
                        vsp_dict = vsp.as_dict()
                        slab = slab_from_file(vsp_path, hkl_tuple)
                        
                        df_list.append(
                            {'hkl': hkl_string, 
                            'hkl_tuple': hkl_tuple,
                            'area': slab.surface_area, 
                            'atoms': vsp_dict['nsites'], 
                            'functional': vsp_dict['run_type'], 
                            'encut': vsp_dict['input']['incar']['ENCUT'], 
                            'ismear': vsp_dict['input']['incar']['ISMEAR'], 
                            'kpoints': vsp_dict['input']['kpoints']['kpoints'],
                            'bandgap': vsp_dict['output']['bandgap'],  
                            'slab_energy': vsp_dict['output']['final_energy'],
                            'slab_per_atom': vsp_dict['output']['final_energy_per_atom'],
                            'vasprun_dictionary': vsp_dict}
                            )

                        if parse_electrostatic: 
                            electrostatic_list.append(
                                parse_electrostatic(path)
                            )
                                                    
                        if parse_core_energy: 
                            core_energy_list.append(
                                parse_core_energy()
                            )
                            
                    
                        
    df = pd.DataFrame(df_list)
    df['surface_energy'] = (
        (df['slab_energy'] - bulk_per_atom * df['atoms'])/(2*df['area']) * 16.02
        ) 

    if electrostatic_list: 
        df['vacuum_potential'] = electrostatic_list
    
    if core_energy_list: 
        df['core_energy'] = core_energy_list

    return df

def parse_electrostatic(path): 
    '''
    checks if a potential.csv file exists in folder, if not tries to calculate
    planar potential from locpot. If neither exist it returns a np.nan 
    '''
    if '{}/potential.csv'.format(path): 
        df = pd.read_csv('{}/potential.csv')
        max_potential = df['planar'].max().round(decimals=3)

    elif '{}/LOCPOT'.format(path): 
        lpt = Locpot.from_file('{}/LOCPOT'.format(path))
        planar = lpt.get_average_along_axis(2)
        max_potential = float(f"{np.max(planar): .3f}")

    else: 
        max_potential = np.nan
        warnings.formatwarning = custom_formatwarning
        warnings.warn('Vacuum electrostatic potential was not parsed - '
        'no LOCPOT or potential.csv files were provided.')

    return max_potential
        

def parse_core_energy(filename=None, core_atom=None, ox_states=None, 
nn_method=CrystalNN()): 
    
    struc = Structure.from_file(filename) 
    struc = oxidation_states(struc)

       
    # need outcar read_core_state_eigen() -> output is dict
    # with [index of atom, zero indexed][str of orbital][ionic run - default=-1 
    # as last run] 
    



                        




            



