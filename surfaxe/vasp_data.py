# Pymatgen  
from pymatgen.core.surface import Slab
from pymatgen import Structure, Specie, Element
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vasp.outputs import Locpot, Outcar, Vasprun
from pymatgen.analysis.local_env import CrystalNN

# Misc
import os
import pandas as pd 
import numpy as np 
import warnings 
warnings.filterwarnings('once')

# surfaxe 
from surfaxe.generation import oxidation_states
from surfaxe.io import _custom_formatwarning, slab_from_file

def process_data(bulk_per_atom, hkl_dict=None, parse_hkl=True, path_to_fols=None,
parse_core_energy=False, core_atom=None, bulk_nn=None, parse_vacuum=True, 
save_csv=True, csv_fname='data.csv', **kwargs): 
    """
    Parses the folders to collect all final data on relevant input and output 
    parameters, and optionally core and vacuum level energies.  

    The function returns None by default and saves the DataFrame to a csv file. 
    Optionally, it can return the DataFrame. 

    Args:
        bulk_per_atom (`float`): Bulk energy per atom in eV per atom. 
        hkl_dict (`dict`, optional): dictionary of tuples of Miller indices 
            and the folders the calculations relating to them are in. Defaults 
            to ``None``. 
            E.g. {(1,-1,2): '1-12'}
        parse_hkl (`bool`, optional): If ``True`` the script parses the names   
            of the folders to get the Miller indices. Defaults to ``True``.
        path_to_fols (`str`, optional): Path to where the hkl folders are. 
            Defaults to None which searches in cwd. 
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
        save_csv (`bool`, optional): If ``True``, it writes data to a csv file.
            Defaults to ``True``.
        csv_fname (`str`, optional): The filename of the csv. Defaults to 
            data.csv 
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
    
    cwd = os.getcwd()
    if path_to_fols: 
        cwd = path_to_fols 

    # Get the Miller indices as tuples and strings from folders in root dir
    if parse_hkl:
        if not hkl_dict: 
            hkl_dict = {}
        for root, fols, files in os.walk(cwd):
            for fol in fols: 
                if fol is not root and len(fol)==3 and fol.isdigit():
                    hkl_dict[tuple(map(int, fol))] = fol

    # Set up additional arguments for get_core_energy 
    get_core_energy_kwargs = {'orbital': '1s', 'ox_states': None, 
    'nn_method': CrystalNN()}
    get_core_energy_kwargs.update(
        (k, kwargs[k]) for k in get_core_energy_kwargs.keys() & kwargs.keys()
    )

    # For each miller index, check if the folders specified are there and 
    # parse them for data
    df_list, electrostatic_list, core_energy_list = ([] for i in range(3))

    for hkl_tuple, hkl_string in hkl_dict.items(): 
        for root, fols, files in os.walk(cwd): 
            for fol in fols: 
                if fol == hkl_string:
                    path = os.path.join(root,fol)
                    vsp_path = '{}/vasprun.xml'.format(path)
                    vsp = Vasprun(vsp_path, parse_potcar_file=False)
                    slab = slab_from_file(vsp.final_structure,hkl_tuple)
                    vsp_dict = vsp.as_dict()
                    
                    df_list.append({
                        'hkl': hkl_string, 
                        'hkl_tuple': hkl_tuple,
                        'area': slab.surface_area, 
                        'atoms': vsp_dict['nsites'], 
                        'functional': vsp_dict['run_type'], 
                        'encut': vsp_dict['input']['incar']['ENCUT'], 
                        'algo': vsp_dict['input']['incar']['ALGO'],
                        'ismear': vsp_dict['input']['parameters']['ISMEAR'],
                        'sigma': vsp_dict['input']['parameters']['SIGMA'],
                        'kpoints': vsp_dict['input']['kpoints']['kpoints'],
                        'bandgap': vsp_dict['output']['bandgap'],  
                        'slab_energy': vsp_dict['output']['final_energy'],
                        'slab_per_atom': vsp_dict['output']['final_energy_per_atom']
                    })

                    if parse_vacuum: 
                        electrostatic_list.append(
                            vacuum(path)
                        )
                                                
                    if parse_core_energy: 
                        otc_path='{}/OUTCAR'.format(path)
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

    # Save to csv or return DataFrame
    if save_csv: 
        df.to_csv(csv_fname, header=True, index=False)
    else:
        return df

def vacuum(path=None): 
    '''
    Gets the energy of the vacuum level. It either parses potential.csv file if 
    available or tries to calculate planar potential from LOCPOT. If neither 
    file is available, function returns np.nan.  

    Args:
        path (`str`, optional): the path to potential.csv or LOCPOT files. 
            Can be the path to a directory in which either file is or you can 
            specify a path that must end in .csv or LOCPOT. Defaults to looking 
            for potential.csv or LOCPOT in cwd. 

    Returns:
        Maximum value of planar potential

    '''
    
    if type(path)==str and path.endswith('.csv'): 
        df = pd.read_csv(path)
        max_potential = df['planar'].max()
        max_potential = round(max_potential, 3)
    
    elif type(path)==str and path.endswith('LOCPOT'): 
        lpt = Locpot.from_file(path)
        planar = lpt.get_average_along_axis(2)
        max_potential = float(f"{np.max(planar): .3f}")
    
    else: 
        if path is None: 
            cwd = os.getcwd()
        else: 
            cwd = path 

        if os.path.isfile('{}/potential.csv'.format(cwd)): 
            df = pd.read_csv('{}/potential.csv'.format(cwd))
            max_potential = df['planar'].max()
            max_potential = round(max_potential, 3)

        elif os.path.isfile('{}/LOCPOT'.format(cwd)): 
            lpt = Locpot.from_file('{}/LOCPOT'.format(cwd))
            planar = lpt.get_average_along_axis(2)
            max_potential = float(f"{np.max(planar): .3f}")

        else: 
            max_potential = np.nan
            warnings.formatwarning = _custom_formatwarning
            warnings.warn('Vacuum electrostatic potential was not parsed - '
            'no LOCPOT or potential.csv files were provided.')

    return max_potential
        

def core_energy(core_atom, bulk_nn, orbital='1s', ox_states=None, 
nn_method=CrystalNN(), outcar='OUTCAR', structure='POSCAR'): 
    """
    Parses the structure and OUTCAR files for the core level energy. Check the 
    validity of nearest neighbour method on the bulk structure before using it 
    on slabs.

    Args: 
        core_atom (`str`, optional): The symbol of atom the core state energy 
            level should be parsed from.  
        bulk_nn (`list`, optional): The symbols of the nearest neighbours of the 
            `core_atom`.   
        orbital (`str`, optional): The orbital of core state. Defaults to 1s.
        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess.  

            Defaults to ``None``.
        nn_method (`class` instance, optional): The coordination number 
            algorithm used. Because the ``nn_method`` is a class, the class 
            needs to be imported from pymatgen.analysis.local_env before it 
            can be instantiated here. Defaults to ``CrystalNN()``.
        outcar (`str`, optional): Path to the OUTCAR file. Defaults to
             ``./OUTCAR``. 
        structure (`str`, optional): Path to the structure file in any format 
            supported by pymatgen. Defaults to ``./POSCAR``. Can also accept a
            pymaten.core.Structure object directly.  

    Returns: 
        Core state energy 
    """
    if type(structure) == str:
        struc = Structure.from_file(structure)
    else:
        struc = structure
    struc = oxidation_states(struc, ox_states)
    bonded_struc = nn_method.get_bonded_structure(struc)
    bulk_nn.sort()
    bulk_nn_str = ' '.join(bulk_nn)
    
    # Get the nearest neighbours info, the c-coordinate and index number 
    # for each atom
    list_of_dicts = [] 
    for n, pos in enumerate(struc): 
        if pos.specie.symbol == core_atom: 
            nn_info = bonded_struc.get_connected_sites(n)
            slab_nn_list = []
            for d in nn_info: 
                nn = d.site.specie.symbol
                slab_nn_list.append(nn)
            slab_nn_list.sort()
            slab_nn = ' '.join(slab_nn_list)
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

    # Check if the df from query isn't empty - if it is it returns
    # a nan as core energy, otherwise it attempts to extract from OUTCAR 
    if type(atom) is np.float64: 
        core_energy = np.nan 
    else:       
        # Read OUTCAR, get the core state energy 
        otc = Outcar(outcar)
        core_energy_dict = otc.read_core_state_eigen()
        try: 
            core_energy = core_energy_dict[atom][orbital][-1]
        except IndexError: 
            core_energy = np.nan



    return core_energy

    



                        




            



