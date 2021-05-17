# Pymatgen
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.analysis.local_env import CrystalNN

# Misc
import os
import pandas as pd
import warnings
import functools
import itertools
import multiprocessing
import copy

# surfaxe
from surfaxe.io import plot_enatom, plot_surfen, slab_from_file,\
 _custom_formatwarning
from surfaxe.vasp_data import vacuum, core_energy

def parse_fols(hkl, bulk_per_atom, path_to_fols=None, parse_core_energy=False,
core_atom=None, bulk_nn=None, parse_vacuum=False, plt_enatom=True, 
plt_enatom_fname='energy_per_atom.png', plt_surfen=True, 
plt_surfen_fname='surface_energy.png', save_csv=True, csv_fname=None, 
verbose=False, **kwargs):
    """
    Parses the convergence folders to get the surface energy, total energy,
    energy per atom, band gap and time taken for each slab and vacuum thickness
    combination. It can optionally parse vacuum and core level energies.

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
            maximum value of planar potential as the vacuum energy level. 
            Defaults to ``False``. 
        plt_enatom (`bool`, optional): Plots the energy per atom. Defaults to 
            ``True``.
        plt_enatom_fname (`str`, optional): The name of the energy per atom plot. 
            Defaults to ``energy_per_atom.png``.
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
    get_core=False
    if parse_core_energy: 
        if core_atom is not None and bulk_nn is not None:
            get_core=True 
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
                functools.partial(_mp_helper, parse_vacuum, parse_core_energy, 
                hkl, core_atom=core_atom, bulk_nn=bulk_nn, 
                **get_core_energy_kwargs), list_of_paths)

        # len(mp_list) == len(list_of_paths), mp_list[0][0] the is main data
        # collected for the dataframe, mp_list[0][1] are the potentials, 
        # mp_list[0][2] are the core energies
        df_list = list(itertools.chain.from_iterable([i[0] for i in mp_list]))
        electrostatic_list = list(itertools.chain.from_iterable(
            [i[1] for i in mp_list]))
        core_energy_list = list(itertools.chain.from_iterable(
            [i[2] for i in mp_list]))
    
    else: 
        df_list, electrostatic_list, core_energy_list = ([] for i in range(3))
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
                    electrostatic_list.append(
                        vacuum(path)
                    )
                                            
            if get_core: 
                core_energy_list.append(
                    core_energy(core_atom, bulk_nn, outcar=otc_path, 
                    structure=slab, **get_core_energy_kwargs)
                    ) 

    df = pd.DataFrame(df_list)
    df['surface_energy'] = (
        (df['slab_energy'] - bulk_per_atom * df['atoms'])/(2*df['area']) * 16.02
        ) 
    
    df['surface_energy_ev'] = (
        (df['slab_energy'] - bulk_per_atom * df['atoms'])/(2*df['area'])
    )

    if electrostatic_list: 
        df['vacuum_potential'] = electrostatic_list
    
    if core_energy_list: 
        df['core_energy'] = core_energy_list

    # Plot energy per atom and surface energy
    plt_kwargs = {'time_taken': True, 'cmap': 'Wistia', 'dpi': 300, 
    'heatmap': False, 'colors': None, 'width': 6, 'height': 5, 'joules': True}
    plt_kwargs.update((k, kwargs[k]) for k in plt_kwargs.keys() & kwargs.keys())

    if plt_surfen: 
        plot_surfen(df, plt_fname=plt_surfen_fname, **plt_kwargs)
    
    if plt_enatom: 
        # Redefine kwargs and delete joules kwarg from the copy
        enatom_kwargs = copy.deepcopy(plt_kwargs) 
        del enatom_kwargs['joules']
        plot_enatom(df, plt_fname=plt_enatom_fname, **enatom_kwargs)

    # Save the csv or return the dataframe
    if save_csv:
        csv_fname = '{}_data.csv'.format(
            ''.join(map(str, hkl))) if csv_fname is None else csv_fname 
        if not csv_fname.endswith('.csv'):
            csv_fname += '.csv'

        df.to_csv(csv_fname, header=True, index=False)
    
    else: 
        return df


def _mp_helper(parse_vacuum, parse_core_energy, hkl, path, slab_thickness,
vac_thickness, slab_index, core_atom=None, bulk_nn=None, **kwargs): 
    """
    Helper function for multiprocessing, returns a list of lists of the main 
    extracted data, electrostatic potential and core energies
    Same args as for parse_fols, only that path is the path to the folder in 
    which the vasprun and OUTCAR for the specific slab/vacuum/index slab are. 
    """
    df_list, electrostatic_list, core_energy_list = ([] for i in range(3))

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
            electrostatic_list.append(
                vacuum(path)
            )
                                    
    if parse_core_energy: 
        core_energy_list.append(
            core_energy(core_atom, bulk_nn, outcar=otc_path, 
            structure=slab, **kwargs)
            ) 

    return [df_list, electrostatic_list, core_energy_list]