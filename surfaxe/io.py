# pymatgen
from pymatgen.io.vasp.sets import DictSet
from pymatgen.core import Structure
from pymatgen.core.surface import Slab

# Misc 
import pandas as pd 
import numpy as np 
import os
import warnings 
import json
from ruamel.yaml import YAML
from pathlib import Path

# Monkeypatching for warnings
def _custom_formatwarning(message, category, filename, lineno, line=''):
    # Ignore everything except the message
    return 'UserWarning: ' + str(message) + '\n'

# Matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler 

   
def slab_from_file(structure, hkl):
    """
    Reads in structure from the file and returns slab object.

    Args:
         structure (str): Structure file in any format supported by pymatgen.
            Will accept a pymatgen.Structure object directly. 
         hkl (tuple): Miller index of the slab in the input file.

    Returns:
         Slab object
    """
    if type(structure) == str:
        slab_input = Structure.from_file(structure)
    else:
        slab_input = structure
    return Slab(slab_input.lattice,
                slab_input.species_and_occu,
                slab_input.frac_coords,
                hkl,
                Structure.from_sites(slab_input, to_unit_cell=True),
                shift=0,
                scale_factor=np.eye(3, dtype=np.int),
                site_properties=slab_input.site_properties)

def slabs_to_file(list_of_slabs, structure, make_fols, make_input_files, 
config_dict, fmt, name, **save_slabs_kwargs): 
    """
    Saves the slabs to file, optionally creates input files. The function can 
    take any relevant keyword argument for DictSet. 

    Args: 
        list_of_slabs (`list`): a list of slab dictionaries made with either of
            surfaxe.generation get_slab functions 
        structure (`str`): Filename of bulk structure file in any format 
            supported by pymatgen. 
        make_fols (`bool`): Makes folders for each termination and slab/vacuum 
            thickness combinations containing structure files. 
            
            * ``True``: A Miller index folder is created, in which folders 
              named slab_vac_index are created to which the relevant structure  
              files are saved. 
                    
                    E.g. for a (0,0,1) slab of index 1 with a slab thickness of 
                    20 Å and vacuum thickness of 30 Å the folder structure would 
                    be: ``001/20_30_1/POSCAR``  

            * ``False``: The indexed structure files are put in a folder named 
              after the bulk formula. 
              
                    E.g. for a (0,0,1) MgO slab of index 1 with a slab thickness 
                    of 20 Å and vacuum thickness of 30 Å the folder structure 
                    would be: ``MgO/POSCAR_001_20_30_1.vasp``
 
        make_input_files (`bool`): Makes INCAR, POTCAR and KPOINTS files in each 
            folder. If ``make_input_files`` is ``True`` but ``make_files`` or 
            ``save_slabs`` is ``False``, files will be saved to folders 
            regardless. This only works with VASP input files, 
            other formats are not yet supported. Defaults to ``False``. 
        config_dict (`dict` or `str`): Specifies the dictionary used for the 
            generation of the input files.
        fmt (`str`, optional): The format of the output files. Options include 
            'cif', 'poscar', 'cssr', 'json', not case sensitive. 
            Defaults to 'poscar'. 
        name (`str`, optional): The name of the surface slab structure file 
            created. Case sensitive. Defaults to 'POSCAR'

    Returns: 
        None, saves surface slabs to file
    """
    bulk_name = structure.composition.reduced_formula

    if make_fols or make_input_files: 
        for slab in list_of_slabs:
            os.makedirs(os.path.join(os.getcwd(), r'{}/{}/{}_{}_{}'.format(
                bulk_name, slab['hkl'], slab['slab_thickness'], 
                slab['vac_thickness'], slab['slab_index'])), exist_ok=True)

            # Makes all input files (KPOINTS, POTCAR, INCAR) based on the config
            # dictionary
            if make_input_files:
                # soft check if potcar directory is set 
                potcars = _check_psp_dir()
                if potcars:
                    cd = _load_config_dict(config_dict)
                    vis = DictSet(slab['slab'], cd, **save_slabs_kwargs)
                    vis.write_input(
                        r'{}/{}/{}_{}_{}'.format(bulk_name, slab['hkl'], 
                        slab['slab_t'], slab['vac_t'],slab['s_index'])
                        )
                # only make the folders with structure files in them
                else: 
                    slab['slab'].to(fmt=fmt,
                filename=r'{}/{}/{}_{}_{}/{}'.format(bulk_name, slab['hkl'],
                slab['slab_t'], slab['vac_t'], slab['s_index'], name))
                warnings.formatwarning = _custom_formatwarning
                warnings.warn('POTCAR directory not set up in pymatgen, only ' 
                'POSCARs were generated ')

            # Just makes the folders with structure files in them
            else:
                slab['slab'].to(fmt=fmt,
                filename=r'{}/{}/{}_{}_{}/{}'.format(bulk_name, slab['hkl'],
                slab['slab_thickness'], slab['vac_thickness'], 
                slab['slab_index'], name))

    # Makes name_hkl_slab_vac_index files in the bulk_name folder
    else:
        suffix='vasp'
        if fmt.lower() != 'poscar': 
            suffix = fmt.lower()
        os.makedirs(os.path.join(os.getcwd(), r'{}'.format(bulk_name)),
        exist_ok=True)
        for slab in list_of_slabs:
            slab['slab'].to(fmt=fmt,
            filename=r'{}/{}_{}_{}_{}_{}.{}'.format(bulk_name, name, 
            slab['hkl'], slab['slab_thickness'], slab['vac_thickness'], 
            slab['slab_index'], suffix))

def _load_config_dict(config_dict=None, path=None): 
    """
    Loads the config dictionary for writing VASP input files. 

    Args: 
        config_dict (``None``, `dict` or `str`): The config dict containing info  
            on INCAR, POTCAR and KPOINTS settings. Can be supplied as: 

            * ``dict``: All settings for the calculations provided as a 
              dictionary of dictionaries 

                    e.g. {'INCAR': {'ENCUT': 500, 'ISYM': 2, 'GGA': 'PE'}, 
                    'KPOINTS': {'reciprocal_density': 20}, 
                    'POTCAR': {'Sn': 'Sn_d', 'O': 'O'}}

            * ``str``: Filename of the config dictionary. Can be a json or a
              yaml file or one of the surfaxe-supplied dictionaries from the
              ``_config_dictionaries`` folder.        

            * ``None``: The default option, makes a PBEsol config dictionary for
              a single shot calculation from the ``PBEsol.json`` file. 
        
        path (``None`` or `str`): The path to where ``_config_dictionaries` are. 
            Needed for CI implementation. 

    Returns: 
        Dictionary
    """

    if path is not None: 
        path_to_conf = path
    else:     
        path_to_conf = str(Path(__file__).parent.joinpath('_config_dictionaries'))

    if type(config_dict) == dict: 
        cd = config_dict

    elif config_dict is not None and type(config_dict)==str: 
        if config_dict.endswith('.json'): 
            with open(config_dict, 'r') as f:
                cd = json.load(f)
        
        elif config_dict.endswith('.yaml'): 
            with open(config_dict, 'r') as y:
                yaml = YAML(typ='safe', pure=True)
                cd = yaml.load(y)
        
        elif 'pe_relax' in config_dict.lower(): 
            with open(os.path.join(path_to_conf, 'PBE_relax.json'), 'r') as f: 
                cd = json.load(f)
        
        elif 'ps_relax' in config_dict.lower(): 
            with open(os.path.join(path_to_conf, 'PBEsol_relax.json'), 'r') as f: 
                cd = json.load(f) 
        
        elif 'pe' in config_dict.lower(): 
            with open(os.path.join(path_to_conf, 'PBE.json'), 'r') as f: 
                cd = json.load(f)
        
        elif 'ps' in config_dict.lower(): 
            with open(os.path.join(path_to_conf, 'PBEsol.json'), 'r') as f: 
                cd = json.load(f)

        elif 'hse' in config_dict.lower(): 
            with open(os.path.join(path_to_conf, 'HSE06.json'), 'r') as f: 
                cd = json.load(f)
        
        else: 
            with open(os.path.join(path_to_conf, 'PBEsol.json'), 'r') as f: 
                cd = json.load(f)
            
            warnings.formatwarning = _custom_formatwarning
            warnings.warn('No valid config dict supplied, reverting to '
            'surfaxe PBEsol.json config dict')
    
    else: 
        with open(os.path.join(path_to_conf, 'PBEsol.json'), 'r') as f: 
            cd = json.load(f)
            
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('No config dict supplied, reverting to surfaxe '
        'PBEsol.json config dict')

    return cd 

def _check_psp_dir(): 
    """
    Helper function to check if potcars are set up correctly for use with 
    pymatgen
    """
    potcar = False
    try: 
        import pymatgen.settings 
        if 'PMG_VASP_PSP_DIR' in pymatgen.settings.SETTINGS:
            potcar = True 
    except ModuleNotFoundError:
        try: 
            import pymatgen 
            if 'PMG_VASP_PSP_DIR' in pymatgen.SETTINGS: 
                potcar = True
        except AttributeError: 
            from pymatgen.core import SETTINGS 
            if 'PMG_VASP_PSP_DIR' in SETTINGS: 
                potcar = True
    return potcar

def _instantiate_structure(structure): 
    """Helper function for instatiating structure files correctly """
    if type(structure) == str:
        struc = Structure.from_file(structure)
    elif type(structure) == Structure or type(structure) == Slab: 
        struc = structure
    else: 
        raise TypeError('structure should either be a file or pmg object')
    
    return struc

def plot_bond_analysis(bond, df=None, filename=None, width=6, height=5, dpi=300,
color=None, plt_fname='bond_analysis.png'): 
    """
    Plots the bond distance with respect to fractional coordinate. Used in 
    conjunction with surfaxe.analysis.bond_analysis.   

    Args:
        bond (`list`): Bond to analyse; e.g. ``['Y', 'O']`` order of elements 
            in the bond must be the same as in the Dataframe or provided file.
        df (`pandas DataFrame`, optional): DataFrame from 
            surfaxe.analysis.bond_analysis. Defaults to ``None``.
        filename (`str`, optional): Path to csv file with data from 
            surfaxe.analysis.bond_analysis. Defaults to ``None``. 
            Either df or filename need to be supplied.
        width (`float`, optional): Width of figure in inches. Defaults to ``6``. 
        height (`float`, optional): Height of figure in inches. Defaults to 
            ``5``. 
        dpi (`int`, optional): Dots per inch. Defaults to ``300``. 
        color (`str`, optional): Color of marker. Defaults to ``None`` which 
            defaults to surfaxe base style
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'bond_analysis.png'``. 
         
    Returns:
        None, saves plot to bond_analysis.png
    """ 
    
    if filename is not None: 
        df = pd.read_csv(filename)
    elif df is not None: 
        df = df
    else: 
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Data not supplied')

    if color is None:
        color = '#F95F6E'

    fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))
    x = df['{}_c_coord'.format(bond[0])]
    y = df['{}-{}_bond_distance'.format(bond[0],bond[1])]
    ax.scatter(x, y, marker='x', markersize=8, c=color)
    ax.set_ylabel("Bond distance / Å")
    ax.legend(['{}-{} bond'.format(bond[0], bond[1])])
    plt.xlabel("Fractional coordinate in c")
    fig.savefig(plt_fname, facecolor='w')
    

def plot_electrostatic_potential(df=None, filename=None, dpi=300, width=6, 
height=5, colors=None, plt_fname='potential.png'): 
    """
    Plots the planar and macroscopic electrostatic potential along one 
    direction. Can take either a DataFrame or a potential.csv file as input. 

    Args: 
        df (`pandas DataFrame`, optional): pandas DataFrame from 
            surfaxe.analysis.electrostatic_potential. Defaults to ``None``. 
        filename (`str`, optional): The filename of csv file with potential 
            data. Defaults to ``None``. 
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        width (`float`, optional): Width of figure in inches. Defaults to ``6``. 
        height (`float`, optional): Height of figure in inches. Defaults to 
            ``5``.
        colors (`list`, optional): A list of colours for planar and macroscopic 
            potential plots. Defaults to ``None``, which defaults to surfaxe 
            base style. 
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'potential.png'``. 
    
    Returns: 
        None, saves plot to potential.png
    """
    if df is not None: 
        df = df
    elif filename is not None: 
        df = pd.read_csv(filename)
    else: 
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Data not supplied')

    if colors==None: 
        colors = ['#FE8995', '#E6505F']

    # Plot both planar and macroscopic, save figure
    fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))
    ax.plot(df['planar'], label='Planar', c=colors[0])
    if 'macroscopic' in df.columns:
        ax.plot(df['macroscopic'], label='Macroscopic', c=colors[1])
    ax.axes.xaxis.set_visible(False)
    ax.legend()
    plt.ylabel('Potential / eV')
    fig.savefig(plt_fname, facecolor='w')

def plot_surfen(df, colors=None, dpi=300, width=8, height=8, 
    plt_fname=None): 
    """
    Plots the surface energy for all terminations. Based on surfaxe.convergence 
    parse_energies. 

    Args: 
        df (pandas DataFrame): DataFrame from `parse_fols`, or any other 
            Dataframe with headings 'slab_thickness', 'vac_thickness', 
            'surface_energy','surface_energy_boettger', 'surface_energy_fm', 
            'time_taken', 'index'. 
        colors (`list`, optional): A list of colours for plots of different 
            surface energies. Defaults to ``None``, which defaults to 
            surfaxe base style. 
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        width (`float`, optional): Width of figure in inches. Defaults to ``8``. 
        height (`float`, optional): Height of figure in inches. Defaults to 
            ``8``.
        plt_fname (`str`, optional): The name of the plot. Defaults to ``None``
            which is either ``surface_energy.png`` for one slab index or 
            ``surface_energy_slab_index.png`` for multiple indices.
    """
    # Make the colour cycler with custom colours
    if colors is None:
        custom_cycler = (cycler(color=['#FE8995','#FE7B88','#E6505F',
            '#CC4755','#9A323C','#802D35','#64232A', '#40161A']))
    else: 
        custom_cycler = (cycler(color=colors))
    
    # make separate plots for different indices 
    for group in df.groupby('slab_index'): 
        df_by_index = group[1].sort_values('vac_thickness')

        # make a 2x2 for 4 thicknesses, 2x3 for 6 and 1xnum for all else 
        num_of_vac = len(df_by_index.groupby('vac_thickness'))
        
        if num_of_vac ==1: 
            fig, ax = plt.subplots(1, 1, dpi=dpi, figsize=(width, height))
            for gp in df_by_index.groupby('vac_thickness'): 
                df_by_vac = gp[1].sort_values('slab_thickness')
                x = df_by_vac['slab_thickness']
                ax.set_prop_cycle(custom_cycler)
                ax.set_xticks(x)
                ax.set_xticklabels(x)
                ax.set_xlabel('Slab thickness / Å')
                ax.set_ylabel('Surface energy / J m$^{-2}$')
                ax.plot(x, df_by_vac['surface_energy'], marker='^', label='Conventional')
                ax.plot(x, df_by_vac['surface_energy_fm'], marker='o', label='Fiorentini-Methfessel')
                ax.plot(x, df_by_vac['surface_energy_boettger'], marker='s', label='Boettger')
                ax.legend()

        elif num_of_vac == 4: 
            fig, ax = plt.subplots(2, 2, dpi=dpi, figsize=(width, height), 
            constrained_layout=True)
            list_of_axes = [ax[0,0], ax[0,1], ax[1,0], ax[1,1]]
            for axs, gp in zip(list_of_axes, df_by_index.groupby('vac_thickness')): 
                df_by_vac = gp[1].sort_values('slab_thickness')
                axs.set_prop_cycle(custom_cycler)
                x = df_by_vac['slab_thickness']
                axs.set_xticks(x)
                axs.set_xticklabels(x)
                axs.set_xlabel('Slab thickness / Å')
                axs.set_ylabel('Surface energy / J m$^{-2}$')
                axs.set_title('Vacuum: {} Å'.format(gp[0]))
                axs.plot(x, df_by_vac['surface_energy'], marker='^', label='Conventional')
                axs.plot(x, df_by_vac['surface_energy_fm'], marker='o', label='Fiorentini-Methfessel')
                axs.plot(x, df_by_vac['surface_energy_boettger'], marker='s', label='Boettger')
                axs.legend()

        elif num_of_vac == 6: 
            fig, ax = plt.subplots(2, 3, dpi=dpi, figsize=(width, height), 
            constrained_layout=True)
            list_of_axes = [ax[0,0], ax[0,1], ax[0,2], ax[1,0], ax[1,1], ax[1,2]]
            for axs, gp in zip(list_of_axes, df_by_index.groupby('vac_thickness')): 
                df_by_vac = gp[1].sort_values('slab_thickness')
                axs.set_prop_cycle(custom_cycler)
                x = df_by_vac['slab_thickness']
                axs.set_xticks(x)
                axs.set_xticklabels(x)
                axs.set_xlabel('Slab thickness / Å')
                axs.set_ylabel('Surface energy / J m$^{-2}$')
                axs.set_title('Vacuum: {} Å'.format(gp[0]))
                axs.plot(x, df_by_vac['surface_energy'], marker='^', label='Conventional')
                axs.plot(x, df_by_vac['surface_energy_fm'], marker='o', label='Fiorentini-Methfessel')
                axs.plot(x, df_by_vac['surface_energy_boettger'], marker='s', label='Boettger')
                axs.legend()
        else: 
            fig, ax = plt.subplots(1, num_of_vac, dpi=dpi, 
                figsize=(width, height), constrained_layout=True)
            for i, gp in enumerate(df_by_index.groupby('vac_thickness')): 
                df_by_vac = gp[1].sort_values('slab_thickness')
                x = df_by_vac['slab_thickness']
                ax[i].set_prop_cycle(custom_cycler)
                ax[i].set_xticks(x)
                ax[i].set_xticklabels(x)
                ax[i].set_xlabel('Slab thickness / Å')
                ax[i].set_ylabel('Surface energy / J m$^{-2}$')
                ax[i].set_title('Vacuum: {} Å'.format(gp[0]))
                ax[i].plot(x, df_by_vac['surface_energy'], marker='^', label='Conventional')
                ax[i].plot(x, df_by_vac['surface_energy_fm'], marker='o', label='Fiorentini-Methfessel')
                ax[i].plot(x, df_by_vac['surface_energy_boettger'], marker='s', label='Boettger')
                ax[i].legend()

        if len(df.groupby('slab_index')) == 1 and plt_fname is None: 
            plt_fname = 'surface_energy.png'
        elif plt_fname is None: 
            plt_fname = 'surface_energy_{}.png'.format(group[0])

        fig.savefig(plt_fname, bbox_inches='tight', facecolor='w')

def plot_enatom(df, colors=None, dpi=300, width=6, height=5, 
plt_fname='energy_per_atom.png'):
    """
    Plots the energy per atom for all terminations. Based on surfaxe.convergence 
    parse_energies.

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols`, or any other 
            Dataframe with headings 'slab_thickness, 'vac_thickness', 
            'slab_per_atom', 'time_taken', 'index'. 
        colors (`list`, optional): A list of colours for plots of different  
            vacuum thicknesses. Defaults to ``None``, which defaults to 
            surfaxe base style. 
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        width (`float`, optional): Width of figure in inches. Defaults to ``6``. 
        height (`float`, optional): Height of figure in inches. Defaults to 
            ``5``. 
        plt_fname (`str`, optional): The name of the plot. Defaults to 
            ``energy_per_atom.png``. If name with no format suffix is supplied,  
            the format defaults to png.

    Returns:
        None, saves energy_per_atom.png
    """
    indices, vals, dfs = ([] for i in range(3))

    # Group the values by termination slab index, create df for time and
    # energy values. Converts the energy and time values to np arrays for 
    # plotting
    for group in df.groupby('slab_index'):
        df2 = group[1].pivot('slab_thickness', 'vac_thickness', 'slab_per_atom')
        indices.append(group[0])
        vals.append(df2.to_numpy())
        dfs.append(df2)
    
    # Make the colour cycler with custom colours
    if colors is None:
        custom_cycler = (cycler(color=['#FE8995','#FE7B88','#E6505F',
            '#CC4755','#9A323C','#802D35','#64232A', '#40161A']))
    else: 
        custom_cycler = (cycler(color=colors))
    
    # Plot only the energy per atom for the only termination present 
    if len(indices)==1: 
        fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))
        for val, df in zip(vals, dfs):
            ax.set_prop_cycle(custom_cycler)
            ax.set_xticks(list(range(len(df.index))))
            ax.set_xticklabels(df.index)
            ax.set_xlabel('Slab thickness / Å')
            ax.set_ylabel('Energy per atom / eV')
            ax.plot(val, marker='x', markersize=8)
            ax.legend(df.columns, title='Vacuum / Å')

    # Plot energy per atom and time taken for different terminations 
    else: 
        fig, ax = plt.subplots(nrows=1, ncols=len(indices), dpi=dpi, 
            figsize=(width, height), constrained_layout=True)
        for i, (index, val, df) in enumerate(zip(indices, vals, dfs)):
            ax[i].set_prop_cycle(custom_cycler)
            ax[i].set_xticks(list(range(len(df.index))))
            ax[i].set_xticklabels(df.index)
            ax[i].set_xlabel('Slab thickness / Å')
            ax[i].set_ylabel('Energy per atom / eV')
            ax[i].plot(val, marker='x', markersize=8)
            ax[i].legend(df.columns, title='Vacuum / Å')
            ax[i].set_title('{}'.format(index))
        
    fig.savefig(plt_fname, bbox_inches='tight', facecolor='w')
