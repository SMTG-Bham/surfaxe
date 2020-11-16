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


# surfaxe 
from surfaxe.convergence import slab_from_file
from surfaxe.generation import oxidation_states
from surfaxe.io import _custom_formatwarning

def data_collection(bulk_per_atom=None, folders=None, hkl_dict=None, 
                    parse_folders=True, parse_core_energy=False, bulk_nn=None,
                    parse_electrostatic=True, save_csv=True, 
                    csv_fname='data.csv',
                    **kwargs): 
    """
    
    Args:
        bulk_per_atom (float): bulk energy per atom from a converged bulk calculation
        folders (list): list of folders where the calculations are ? like pbesol
        hse06 - whatever the folders that contain the calcs are called  
        hkl_dict (dict): dictionary of tuples of Miller indices and the folders the 
        calculations relating to them are in, e.g. {(1,-1,2): '1-12'}; default=None
        parse_fols (bool): if true the script parses the root folders to get the 
        Miller indices; default=True.
        parse_core_energy (bool): if True will attempt to parse core energies 
        from a supplied OUTCAR for core level energies; default=False.
        bulk_nn (list, optional): 
        parse_electrostatic (bool): if True will attempt to parse LOCPOT using 
        analysis.electrostatic_potential for vacuum potential; default=True. 
        save_csv (bool): whether or not to save to csv; default=True.
        csv_fname (str): filename of the csv; default='data.csv' 
    Returns: 
        DataFrame
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
        if not hkl_dict: 
            hkl_dict = {}
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
                        
                        df_list.append({
                            'hkl': hkl_string, 
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
                            'vasprun_dictionary': vsp_dict
                        })

                        if parse_electrostatic: 
                            electrostatic_list.append(
                                get_vacuum_level(path)
                            )
                                                    
                        if parse_core_energy: 
                            core_energy_list.append(
                                get_core_energy(path, bulk_nn)
                            )
                            
                    
                        
    df = pd.DataFrame(df_list)
    df['surface_energy'] = (
        (df['slab_energy'] - bulk_per_atom * df['atoms'])/(2*df['area']) * 16.02
        ) 

    if electrostatic_list: 
        df['vacuum_potential'] = electrostatic_list
    
    if core_energy_list: 
        df['core_energy'] = core_energy_list

    # Save to csv or return DataFrame
    if save_csv: 
        save_csv(df, **kwargs)
    else:
        return df

def get_vacuum_level(path): 
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
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Vacuum electrostatic potential was not parsed - '
        'no LOCPOT or potential.csv files were provided.')

    return max_potential
        

def get_core_energy(path, bulk_nn, orbital='1s', ox_states=None, 
nn_method=CrystalNN()): 
    """
    Args: 

    Returns: 
        core state energy 
    """
    struc = Structure.from_file('{}/vasprun.xml'.format(path)) 
    struc = oxidation_states(struc)
    bulk_nn.sort()
    bulk_nn_str = ''.join(bulk_nn)
    
    # Get the nearest neighbours info, the c-coordinate and index number 
    # for each atom
    list_of_dicts = [] 
    for n, pos in enumerate(struc): 
        nn_info = nn_method.get_nn_info(struc, n)
        slab_nn_list = []
        for d in nn_info: 
            nn = d.get('site').specie.symbol
            slab_nn_list.append(nn)
        slab_nn_list.sort()
        slab_nn = ''.join(slab_nn_list)
        list_of_dicts.append({
            'atom': n, 
            'nn': slab_nn, 
            'c_coord': pos.c
        })

    # Make pandas Dataframe, query the interquartile range of fractional 
    # coordinates of atoms whose nearest neighbour environment is same as the
    # bulk nearest neighbours provided, get the atom that is the nearest to the 
    # median of the interquartile range 
    df = pd.DataFrame(list_of_dicts)
    low, high = df['c_coord'].quantile([0.25,0.75])
    df = df.query('@low<c_coord<@high and nn==@bulk_nn_str')
    atom = df['atom'].quantile(interpolation='nearest')        

    # Read OUTCAR, get the core state energy 
    otc = Outcar('{}/OUTCAR'.format(path))
    core_energy_dict = otc.read_core_state_eigen()
    core_energy = core_energy_dict[atom][orbital][-1]

    return core_energy

    



                        




            



