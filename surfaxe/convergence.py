# Pymatgen
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.analysis.local_env import CrystalNN

# Misc
import os
import pandas as pd
import numpy as np
import warnings
import functools
import itertools
import multiprocessing
import copy
from sklearn.linear_model import LinearRegression
import json

# surfaxe
from surfaxe.io import plot_surfen, slab_from_file, _custom_formatwarning
from surfaxe.vasp_data import vacuum, core_energy
from surfaxe.analysis import bond_analysis, electrostatic_potential

# todo: 
# - add a script that takes json from here and generate and does cart disp 
# between the two them 
# - add docs for parse_structures 

def parse_energies(hkl, bulk_per_atom, path_to_fols=None, parse_core_energy=False,
core_atom=None, bulk_nn=None, parse_vacuum=False, remove_first_energy=False,
plt_surfen=True, plt_surfen_fname='surface_energy.png', save_csv=True,
csv_fname=None, verbose=False, **kwargs):
    """
    Parses the convergence folders to get the surface energy, total energy,
    energy per atom, band gap and time taken for each slab and vacuum thickness
    combination. Calculates surface energies using Fiorentini-Methfessel and 
    Boettger methods. It can optionally parse vacuum and core level energies.

    ``path_to_fols`` specifies the parent directory containing subdirectories 
    that must include the miller index specified. e.g. if ``hkl=(0,0,1)`` there 
    must be a ``001/`` subdirectory present somewhere on the path. Each 
    directory within the subdirectory must contain a vasprun.xml and OUTCAR file.  

    Args:
        hkl (`tuple`): Miller index of the slab.
        bulk_per_atom (`float`): Bulk energy per atom from a converged 
            bulk calculation in eV per atom.
        path_to_fols (`str`, optional): Path to the convergence folders. 
            Defaults to None which is cwd
        parse_core_energy (`bool`, optional): If ``True`` the script attempts to 
            parse core energies from a supplied OUTCAR. Defaults to ``False``. 
        core_atom (`str`, optional): The symbol of atom the core state energy 
            level should be parsed from. Defaults to ``None``. 
        bulk_nn (`list`, optional): The symbols of the nearest neighbours of the 
            `core_atom`. Defaults to ``None``. 
        parse_vacuum (`bool`, optional): if ``True`` the script attempts 
            to parse LOCPOT using analysis.electrostatic_potential to use the 
            max value of planar potential as the vacuum energy level (eV). 
            It also calculates the average gradient of the vacuum region (meV). Defaults to ``False``. 
        remove_first_energy (`bool`, optional): Remove the first data point in 
            calculation of Fiorentini-Metfessel and Boettger surface energy. 
            Use if the first energy is somewhat of an outlier.
        plt_surfen (`bool`, optional): Plots the surface energy. Defaults to 
            ``True``.
        plt_surfen_fname (`str`, optional): The name of the surface energy plot.
            Defaults to ``surface_energy.png``
        save_csv (`bool`, optional): Saves the csv. Defaults to ``True``.
        csv_fname (`str`, optional): Name of the csv file to save. Defaults to
            hkl_data.csv, where hkl are the miller indices.
        verbose (`bool`, optional): Whether or not to print extra info about the
            folders being parsed. Defaults to ``False``. 

    Returns:
        DataFrame 
    """
    
    # Update kwargs for core energy 
    get_core_energy_kwargs = {'orbital': '1s', 'ox_states': None, 
    'nn_method': CrystalNN()}
    get_core_energy_kwargs.update(
        (k, kwargs[k]) for k in get_core_energy_kwargs.keys() & kwargs.keys()
    )
    get_core = False
    if parse_core_energy: 
        if core_atom is not None and bulk_nn is not None:
            get_core = True 
        else: 
            warnings.formatwarning = _custom_formatwarning
            warnings.warn(('Core atom or bulk nearest neighbours were not '
            'supplied. Core energy will not be parsed.'))

    # Set directory 
    cwd = os.getcwd() if path_to_fols is None else path_to_fols

    # Get all paths to slab_vac_index folders
    list_of_paths=[]
    for root, fols, files in os.walk(cwd):
        for fol in fols:
            # Perform a loose check that we are looking in the right place, 
            # also avoid .ipynb_checkpoint files 
            if ''.join(map(str,hkl)) in root.split('/') and '.' not in fol:
                if len(fol.split('_')) == 3:
                    list_of_paths.append([
                        os.path.join(root, fol),
                        fol.split('_')[0], 
                        fol.split('_')[1], 
                        fol.split('_')[2]
                    ])
                    if verbose:
                        print(root, fol)
    
    if len(list_of_paths) > 20 and parse_core_energy:
        warnings.formatwarning = _custom_formatwarning
        warnings.warn(('Determining core energies for {} slabs may be slow. ' 
        'Running on {} cores.').format(len(list_of_paths), 
                                       multiprocessing.cpu_count()))


    # Check if multiple cores are available, iterate through paths to folders 
    # and parse folders 
    if multiprocessing.cpu_count() > 1: 
        with multiprocessing.Pool() as pool: 
            mp_list = pool.starmap(
                functools.partial(_mp_helper_energy, parse_vacuum, 
                get_core, hkl, core_atom=core_atom, bulk_nn=bulk_nn, 
                **get_core_energy_kwargs), list_of_paths)

        # len(mp_list) == len(list_of_paths), mp_list[0][0] the is main data
        # collected for the dataframe, mp_list[0][1] are the potentials, 
        # mp_list[0][2] are the core energies
        df_list = list(itertools.chain.from_iterable([i[0] for i in mp_list]))
        electrostatic_list = list(itertools.chain.from_iterable(
            [i[1] for i in mp_list]))
        core_energy_list = list(itertools.chain.from_iterable(
            [i[2] for i in mp_list]))
        gradient_list = list(itertools.chain.from_iterable(
            [i[3] for i in mp_list]))
    
    else: 
        df_list, electrostatic_list, core_energy_list, gradient_list = ([] for i in range(4))
        for path, slab_thickness, vac_thickness, slab_index in list_of_paths: 
            vsp_path = '{}/vasprun.xml'.format(path)
            otc_path = '{}/OUTCAR'.format(path)

            # instantiate structure, slab, vasprun and outcar objects
            vsp = Vasprun(vsp_path, parse_potcar_file=False)
            otc = Outcar(otc_path)
            slab = slab_from_file(vsp.final_structure, hkl)
            vsp_dict = vsp.as_dict()

            # extract the time data
            otc_times = otc.run_stats

            df_list.append(
                {'hkl_string': ''.join(map(str, hkl)), 
                'hkl_tuple': hkl, 
                'slab_thickness': slab_thickness,
                'vac_thickness': vac_thickness,
                'slab_index': slab_index,
                'atoms': vsp_dict['nsites'], 
                'area': slab.surface_area, 
                'bandgap': vsp_dict['output']['bandgap'],
                'slab_energy': vsp_dict['output']['final_energy'],
                'slab_per_atom': vsp_dict['output']['final_energy_per_atom'],
                'time_taken': otc_times['Elapsed time (sec)']})

            if parse_vacuum: 
                # get value of potential in vacuum
                electrostatic_list.append(vacuum(path))

                lpt = '{}/LOCPOT'.format(path) 
                p = electrostatic_potential(lpt, save_plt=False, save_csv=False)
                # gradient in eV 
                g = p['gradient'].to_numpy()
                ratio = int(vac_thickness)/(int(vac_thickness)+int(slab_thickness))
                # check if slab is centred so can get the gradient
                if 0.45 < slab.center_of_mass[2] < 0.55:
                    # divide by 2 bc two regions, 0.75 is a scaling factor that 
                    # should hopefully work to avoid any charge from dangling 
                    # bonds? 
                    # mean is in meV 
                    a = int(len(g)*ratio/2*0.75)
                    gradient_list.append(np.mean(g[:a]) * 1000) 

                else: 
                    a = int(len(g)*ratio*0.75)
                    # if atoms are more towards the end of the slab, the start 
                    # should be vacuum
                    if slab.center_of_mass[2] > 0.5: 
                        gradient_list.append(np.mean(g[:a]) * 1000)  # meV 
                    # and then if atoms are at the start, the vacuum is at
                    # the end
                    else: 
                        gradient_list.append(np.mean(g[(len(g)-a):]) * 1000) 

            if get_core: 
                core_energy_list.append(
                    core_energy(core_atom, bulk_nn, outcar=otc_path, 
                    structure=slab, **get_core_energy_kwargs)
                    ) 

    df = pd.DataFrame(df_list)

    if electrostatic_list: 
        df['vacuum_potential'] = electrostatic_list
    
    if gradient_list: 
        df['vacuum_gradient'] = gradient_list # in meV
    
    if core_energy_list: 
        df['core_energy'] = core_energy_list
        
    df['surface_energy'] = (
        (df['slab_energy'] - bulk_per_atom * df['atoms'])/(2*df['area']) * 16.02
        ) 
    
    # Add Fiorentini-Methfessel and Boettger methods for calculating 
    # surface energies 
    dfs = []
    for index in df.groupby('slab_index'): 
        df_index = index[1].sort_values('vac_thickness')
        for group in df_index.groupby('vac_thickness'): 
            df2 = group[1].sort_values('slab_thickness')
            # deepcopy to keep it fresh 
            df3 = copy.deepcopy(df2)
            # Fiorentini-methfessel
            # remove first data point if it's too much of an outlier and there are 
            # more than three data points 
            if remove_first_energy and len(df2['atoms']) >= 3: 
                x = np.delete(df2['atoms'].to_numpy(), 0).reshape(-1,1)
                y = np.delete(df2['slab_energy'].to_numpy(), 0) 
            else: 
                x = df2['atoms'].to_numpy().reshape(-1,1)
                y = df2['slab_energy'].to_numpy()
            
            model = LinearRegression().fit(x,y)
            df3['surface_energy_fm'] = (
            (df3['slab_energy'] - model.coef_ * df3['atoms'])/(2*df3['area'])*16.02)
            
            # Boettger
            # new dataframe, just in case 
            df4 = group[1].sort_values('slab_thickness')
            if remove_first_energy and len(df4['atoms']) >= 3: 
                # big energy = M+1 layers energy, small energy = M layers, 
                # calculating bulk energy for M layers
                # remove first two from M+1 energy, add a nan in place of one and 
                # at the end; also replace the first from the M energies with nan; 
                # same with atoms
                big_energy = df4['slab_energy'].iloc[2:].to_numpy()
                big_energy = np.append(big_energy, np.nan)
                big_energy = np.insert(big_energy, 0, np.nan)
                
                small_energy = df4['slab_energy'].iloc[1:].to_numpy()
                small_energy = np.insert(small_energy, 0, np.nan)
                
                big_atoms = df4['atoms'].iloc[2:].to_numpy()
                big_atoms = np.append(big_atoms, np.nan)
                big_atoms = np.insert(big_atoms, 0, np.nan)
                small_atoms = df4['atoms'].astype('float').iloc[1:].to_numpy()
                small_atoms = np.insert(small_atoms, 0, np.nan) 
            
            else: 
                # need to remove the first from big energy and add a nan at the end 
                big_energy = df4['slab_energy'].iloc[1:].to_numpy()
                big_energy = np.append(big_energy, np.nan)
                
                small_energy = df4['slab_energy'].to_numpy()
                
                big_atoms = df4['atoms'].iloc[1:].to_numpy()
                big_atoms = np.append(big_atoms, np.nan)
                small_atoms = df4['atoms'].to_numpy()
            
            # difference M+1 - M, get bulk energy from E(M+1)-E(M) / (M+1)-M        
            diff_energy = big_energy-small_energy
            diff_atoms = big_atoms-small_atoms
            bulk_energies = diff_energy/diff_atoms
            # make last energy a nan to get a nan surface energy 

            small_energy[-1] = np.nan
        
            df3['surface_energy_boettger'] = (
            (small_energy - df3['atoms']*bulk_energies) / (2*df3['area']) * 16.02)

            dfs.append(df3)
    
    # Concat list back to one dataframe
    df = pd.concat(dfs)

    if remove_first_energy and len(df2['atoms']) < 3:
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('First data point was not removed - less than three '
        'data points were present in dataset') 

    # Plot surface energy
    plt_kwargs = {'colors': None, 'width': 6, 'height': 5}
    plt_kwargs.update((k, kwargs[k]) for k in plt_kwargs.keys() & kwargs.keys())

    if plt_surfen: 
        plot_surfen(df, plt_fname=plt_surfen_fname, **plt_kwargs)

    # Save the csv or return the dataframe
    if save_csv:
        csv_fname = '{}_data.csv'.format(
            ''.join(map(str, hkl))) if csv_fname is None else csv_fname 
        if not csv_fname.endswith('.csv'):
            csv_fname += '.csv'

        df.to_csv(csv_fname, header=True, index=False)
    
    else: 
        return df

def parse_structures(hkl, structure_file='CONTCAR', bond=None, nn_method=CrystalNN(), 
 path_to_fols=None, verbose=False, json_fname=None,  **kwargs): 
    # collect the relaxed structures & put them in a json file 
    # has to be same format as the input (-layers) 

    # Set directory 
    cwd = os.getcwd() if path_to_fols is None else path_to_fols

    # Get all paths to slab_vac_index folders, list=[[path,slab,vac,index],..]
    list_of_paths=[]
    for root, fols, files in os.walk(cwd):
        for fol in fols:
            # Perform a loose check that we are looking in the right place, 
            # also avoid .ipynb_checkpoint files 
            if ''.join(map(str,hkl)) in root.split('/') and '.' not in fol:
                if len(fol.split('_')) == 3:
                    list_of_paths.append([
                        os.path.join(root, fol),
                        fol.split('_')[0], 
                        fol.split('_')[1], 
                        fol.split('_')[2]
                    ])
                    if verbose:
                        print(root, fol)

    lst = []
    for path, slab_thickness, vac_thickness, slab_index in list_of_paths: 
        struc_path = '{}/{}'.format(path, structure_file)
        slab = slab_from_file(struc_path, hkl)
        lst.append({
            'hkl': hkl, 
            'slab_thickness': slab_thickness, 
            'vac_thickness': vac_thickness, 
            'slab_index': slab_index,
            'slab': slab.as_dict()
        })
        # if bond is set, do bond analysis; single bond
        if bond is not None and type(bond[0]) == str: 
            csv_fname = '{}_bond_analysis_{}.csv'.format(path, ''.join(bond))
            bond_analysis(slab, bond, nn_method=nn_method, 
            csv_fname=csv_fname, **kwargs)
            
        # multiple bonds 
        elif bond is not None and type(bond[0]) == list: 
            for b in bond: 
                csv_fname = '{}_bond_analysis_{}.csv'.format(path, ''.join(b))
                bond_analysis(slab, b, nn_method=nn_method, 
                csv_fname=csv_fname, **kwargs)

    if json_fname is None: 
        bulk_name = slab.composition.reduced_formula
        json_fname = '{}_parsed_metadata.json'.format(bulk_name)
    
    with open(json_fname, 'w') as f: 
        json.dump(lst, f)



def _mp_helper_energy(parse_vacuum, get_core, hkl, path, slab_thickness,
vac_thickness, slab_index, core_atom=None, bulk_nn=None, **kwargs): 
    """
    Helper function for multiprocessing, returns a list of lists of the main 
    extracted data, electrostatic potential and core energies
    Same args as for parse_fols, only that path is the path to the folder in 
    which the vasprun and OUTCAR for the specific slab/vacuum/index slab are. 
    """
    df_list, electrostatic_list, core_energy_list, gradient_list = ([] for i in range(4))

    # instantiate structure, slab, vasprun and outcar objects
    vsp_path = '{}/vasprun.xml'.format(path)
    otc_path = '{}/OUTCAR'.format(path)
    vsp = Vasprun(vsp_path, parse_potcar_file=False)
    otc = Outcar(otc_path)
    slab = slab_from_file(vsp.final_structure, hkl)
    vsp_dict = vsp.as_dict()

    # extract the time data
    otc_times = otc.run_stats

    df_list.append(
        {'hkl_string': ''.join(map(str, hkl)), 
        'hkl_tuple': hkl, 
        'slab_thickness': slab_thickness,
        'vac_thickness': vac_thickness,
        'slab_index': slab_index,
        'atoms': vsp_dict['nsites'], 
        'area': slab.surface_area, 
        'bandgap': vsp_dict['output']['bandgap'],
        'slab_energy': vsp_dict['output']['final_energy'],
        'slab_per_atom': vsp_dict['output']['final_energy_per_atom'],
        'time_taken': otc_times['Elapsed time (sec)']})

    if parse_vacuum: 
        electrostatic_list.append(vacuum(path))

        lpt = '{}/LOCPOT'.format(path) 
        p = electrostatic_potential(lpt, save_plt=False, save_csv=False)
        # gradient in eV 
        g = p['gradient'].to_numpy()
        ratio = int(vac_thickness)/(int(vac_thickness)+int(slab_thickness))
        # check if slab is centred so can get the gradient
        if 0.45 < slab.center_of_mass[2] < 0.55:
            # divide by 2 bc two regions, 0.75 is a scaling factor that 
            # should hopefully work to avoid any charge from dangling 
            # bonds? 
            a = int(len(g)*ratio/2*0.75)
            gradient_list.append(np.mean(g[:a])* 1000) 

        else: 
            a = int(len(g)*ratio*0.75)
            # if atoms are more towards the end of the slab, the start 
            # should be vacuum
            if slab.center_of_mass[2] > 0.5: 
                gradient_list.append(np.mean(g[:a])* 1000)  
            # and then if atoms are at the start, the vacuum is at
            # the end
            else: 
                gradient_list.append(np.mean(g[(len(g)-a):]) * 1000)  
                                    
    if get_core: 
        core_energy_list.append(
            core_energy(core_atom, bulk_nn, outcar=otc_path, 
            structure=slab, **kwargs)
            ) 

    return [df_list, electrostatic_list, core_energy_list, gradient_list]