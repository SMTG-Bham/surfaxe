# pymatgen
from pymatgen.io.vasp.sets import DictSet
from pymatgen import Structure
from pymatgen.core.surface import Slab

# Misc 
import pandas as pd 
import numpy as np 
import os
import warnings 
import json
from pathlib import Path

# Monkeypatching for warnings
def _custom_formatwarning(message, category, filename, lineno, line=''):
    # Ignore everything except the message
    return 'UserWarning: ' + str(message) + '\n'

# Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler 

def load_config_dict(config_dict, path_to_config_dir=None): 
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

            * ``str``: Filename of the config dictionary in the 
              ``_config_directories`` folder. If the filename does not exist,  
              the function defaults to the ``PBEsol_config.json`` file.            

            * ``None``: The default option, makes a PBEsol config dictionary for
              a single shot calculation from the ``PBEsol_config.json`` file. 
        path_to_config_dir (`str`, optional). The path to the directory in which 
            the json config files are. Defaults to 
            ``surfaxe/surfaxe/_config_dictionaries``. 
    Returns: 
        Dictionary
    
    """

    
    if type(path_to_config_dir) is str: 
        conf_dir = path_to_config_dir
    else: 
        conf_dir = str(Path(__file__).parent.joinpath('_config_dictionaries'))

    if type(config_dict) is dict: 
        cd = config_dict 
    elif type(config_dict) is str: 
        if os.path.isfile(os.path.join(conf_dir, config_dict)): 
            with open(os.path.join(conf_dir, config_dict), 'r') as f:
                cd = json.load(f)
        else: 
            with open(os.path.join(conf_dir, 'PBEsol_config.json'), 'r') as f:
                cd = json.load(f)
    else: 
        with open(os.path.join(conf_dir, 'PBEsol_config.json'), 'r') as f:
            cd = json.load(f)

    return cd 

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
    struc = Structure.from_file(structure)
    bulk_name = struc.formula.replace(" ", "")

    if make_fols or make_input_files: 
        for slab in list_of_slabs:
            os.makedirs(os.path.join(os.getcwd(), r'{}/{}/{}_{}_{}'.format(
                bulk_name, slab['hkl'],
            slab['slab_t'], slab['vac_t'], slab['s_index'])), exist_ok=True)

            # Makes all input files (KPOINTS, POTCAR, INCAR) based on the config
            # dictionary
            if make_input_files:
                cd = load_config_dict(config_dict)
                vis = DictSet(slab['slab'], cd, **save_slabs_kwargs)
                vis.write_input(
                    r'{}/{}/{}_{}_{}'.format(bulk_name, slab['hkl'], 
                    slab['slab_t'], slab['vac_t'],slab['s_index'])
                    )

            # Just makes the folders with structure files in them
            else:
                slab['slab'].to(fmt=fmt,
                filename=r'{}/{}/{}_{}_{}/{}'.format(bulk_name, slab['hkl'],
                slab['slab_t'], slab['vac_t'], slab['s_index'], name))

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
            slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index'], suffix))

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
            ``'bond_analysis.png'``. If name with no format suffix is supplied,  
            the format defaults to png.
         
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
        color = '#FF596A'
    
    if not plt_fname.endswith('.png'):
        plt_fname += '.png'

    fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))
    x = df['{}_c_coord'.format(bond[0])]
    y = df['{}-{}_bond_distance'.format(bond[0],bond[1])]
    ax.scatter(x, y, marker='x', c=color)
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
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'potential.png'``. If name with no format suffix is supplied,  
            the format defaults to png.
        colors (`list`, optional): A list of colours for planar and macroscopic 
            potential plots. Defaults to ``None``, which defaults to surfaxe 
            base style. 
    
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

    if colors==None or len(colors)<2: 
        colors = ['#FFADB6', '#FF596A']
    
    if not plt_fname.endswith('.png'):
        plt_fname += '.png'

    # Plot both planar and macroscopic, save figure
    fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))
    ax.plot(df['planar'], label='Planar', c=colors[0])
    ax.plot(df['macroscopic'], label='Macroscopic', c=colors[1])
    ax.axes.xaxis.set_visible(False)
    ax.legend()
    plt.ylabel('Potential / eV')
    fig.savefig(plt_fname, facecolor='w')

def plot_surfen(df, time_taken=True, colors=None, dpi=300, width=6, height=5, 
heatmap=False, cmap='Wistia',  plt_fname='surface_energy.png'):
    """
    Plots the surface energy for all terminations. Based on surfaxe.convergence 
    parse_fols. 

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols`, or any other 
            Dataframe with headings 'slab_thickness', 'vac_thickness', 
            'surface_energy', 'time_taken', 'index'. 
        time_taken (bool): Show the time taken for calculation to finish on the
            figure. Defaults to True.
        colors (`list`, optional): A list of colours for plots of different 
            vacuum thicknesses. Defaults to ``None``, which defaults to 
            surfaxe base style. 
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        width (`float`, optional): Width of figure in inches. Defaults to ``6``. 
        height (`float`, optional): Height of figure in inches. Defaults to 
            ``5``. 
        heatmap (`bool`, optional): If True plots a heatmap of surface energies.
            Defaults to False.
        cmap (`str`, optional): Matplotlib colourmap. Defaults to 'Wistia'
        plt_fname (`str`, optional): The name of the plot. Defaults to 
            ``surface_energy.png``. If name with no format suffix is supplied,  
            the format defaults to png.
        
    Returns:
        None, saves surface_energy.png to file
    """
    if not plt_fname.endswith('.png'):
        plt_fname += '.png'
    
    indices, vals, times, dfs, dfs_times = ([] for i in range(5))

    # Group the values by termination slab index, create df for time and
    # energy values. Converts the energy and time values to np arrays for 
    # plotting
    for group in df.groupby('slab_index'):
        df2 = group[1].pivot('slab_thickness', 'vac_thickness', 'surface_energy')
        df3 = group[1].pivot('slab_thickness', 'vac_thickness', 'time_taken')
        indices.append(group[0])
        vals.append(df2.to_numpy())
        times.append(df3.to_numpy())
        dfs.append(df2)
        dfs_times.append(df3)

    if heatmap:
        # Plotting has to be separated into plotting for one and more than one
        # indices mpl won't let you index axes if there's only one set
        if len(indices) == 1:
            fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))

            for (index, val, time, df) in zip(indices, vals, times, dfs):
                ax.set_yticks(list(range(len(df.index))))
                ax.set_yticklabels(df.index)
                ax.set_ylabel('Slab thickness / Å')
                ax.set_xticks(list(range(len(df.columns))))
                ax.set_xticklabels(df.columns)
                ax.set_xlabel('Vacuum thickness / Å')
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.2)
                im = ax.imshow(val, cmap=cmap)
                cbar = fig.colorbar(im, cax=cax, orientation='vertical')
                cbar.set_label('Surface energy / J m$^{-2}$')
                ax.invert_yaxis()

            # Add the surface energy value labels to the plot - the for loops 
            # are needed in this order and they can't be zipped because this 
            # is the only way they display the values correctly
            for df in dfs:
                for j in range(len(df.index)):
                    for k in range(len(df.columns)):
                        for val in vals:
                            ax.text(k, j, f"{val[j, k]: .3f}", ha="center",
                                            va="bottom", color="black")

            # Add the time taken labels to the plot
            if time_taken:
                for df in dfs:
                    for j in range(len(df.index)):
                        for k in range(len(df.columns)):
                            for time in times:
                                ax.text(k, j, (f"{time[j, k]: .0f}"+' s'),
                                ha="center", va="top", color="black")

        # Plotting for multiple indices
        else:
            fig, ax = plt.subplots(ncols=len(indices), dpi=dpi, 
            figsize=(width, height))

            # Iterate through the values for plotting, create each plot on a 
            # separate ax, add the colourbar to each ax
            for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):
                ax[i].set_title('{}'.format(index))
                ax[i].set_yticks(list(range(len(df.index))))
                ax[i].set_yticklabels(df.index)
                ax[i].set_ylabel('Slab thickness / Å')
                ax[i].set_xticks(list(range(len(df.columns))))
                ax[i].set_xticklabels(df.columns)
                ax[i].set_xlabel('Vacuum thickness / Å')
                im = ax[i].imshow(val, cmap=cmap)
                divider = make_axes_locatable(ax[i])
                cax = divider.append_axes("right", size="5%", pad=0.2)
                cbar = plt.colorbar(im, cax=cax)
                cbar.set_label('Surface energy / J m$^{-2}$')
                ax[i].invert_yaxis()
            fig.tight_layout()

            # Add the surface energy value labels to the plot
            for df in dfs:
                for j in range(len(df.index)):
                    for k in range(len(df.columns)):
                        for i, val in enumerate(vals):
                            ax[i].text(k, j, f"{val[j, k]: .3f}", ha="center",
                                            va="bottom", color="black")

            # Add the time taken labels to the plot
            if time_taken:
                for df in dfs:
                    for j in range(len(df.index)):
                        for k in range(len(df.columns)):
                            for i, time in enumerate(times):
                                ax[i].text(k, j, (f"{time[j, k]: .0f}"+' s'),
                                ha="center", va="top", color="black")
    # Line plots
    else: 
        # Make the colour cycler with custom colours
        if colors is None:
            custom_cycler = (cycler(color=['#FFADB6','#FF596A','#CC4755',
                '#802D35','#E6505F','#40161A','#CC858C']))
        else: 
            custom_cycler = (cycler(color=colors))

        # Get the number of subplots 
        nrows = len(indices)
        ncols = 1
        if time_taken: 
            ncols = 2
        
        # Plot only the surface energy for the only termination present 
        if nrows==1 and ncols==1: 
            fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))
            for index, val, time, df in zip(indices, vals, times, dfs):
                ax.set_prop_cycle(custom_cycler)
                ax.set_title(index)
                ax.set_xticks(list(range(len(df.index))))
                ax.set_xticklabels(df.index)
                ax.set_xlabel('Slab thickness / Å')
                ax.set_ylabel('Surface energy / J m$^{-2}$')
                ax.plot(val, marker='x')
                ax.legend(df.columns, title='Vacuum / Å')

        # Plot surface energy and time taken for the only termination present
        elif nrows==1 and ncols==2: 
            fig, ax = plt.subplots(1, 2, dpi=dpi, figsize=(width, height))
            for index, val, time, df in zip(indices, vals, times, dfs):
                ax[0].set_prop_cycle(custom_cycler)
                ax[0].set_title(index)
                ax[0].set_xticks(list(range(len(df.index))))
                ax[0].set_xticklabels(df.index)
                ax[0].set_xlabel('Slab thickness / Å')
                ax[0].set_ylabel('Surface energy / J m$^{-2}$')
                ax[0].plot(val, marker='x')
                ax[0].legend(df.columns, title='Vacuum / Å')
                
                ax[1].set_prop_cycle(custom_cycler)
                ax[1].set_title(index)
                ax[1].set_xticks(list(range(len(df.index))))
                ax[1].set_xticklabels(df.index)
                ax[1].set_xlabel('Slab thickness / Å')
                ax[1].set_ylabel('Time taken / s')
                ax[1].plot(time, marker='x')
                ax[1].legend(df.columns, title='Vacuum / Å')
                plt.tight_layout()

        # Plot surface energies and or times taken for different terminations 
        else: 
            fig, ax = plt.subplots(nrows=nrows,ncols=ncols, dpi=dpi, 
            figsize=(width, height))
            for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):    

                # Separate plotting of times and energies, energies in the first
                # column, times in the second. Annoyingly matplotlib doesn't 
                # like ax[i,0] if there isn't a second axis so need to separate 
                # it for time_taken True and False
            
                if time_taken: 
                    ax[i,0].set_prop_cycle(custom_cycler)
                    ax[i,0].set_xticks(list(range(len(df.index))))
                    ax[i,0].set_xticklabels(df.index)
                    ax[i,0].set_xlabel('Slab thickness / Å')
                    ax[i,0].set_ylabel('Surface energy / J m$^{-2}$')
                    ax[i,0].plot(val, marker='x')
                    ax[i,0].legend(df.columns, title='Vacuum / Å')
                    ax[i,0].set_title('{}'.format(indices[i]))

                    ax[i,1].set_prop_cycle(custom_cycler)
                    ax[i,1].set_xticks(list(range(len(df.index))))
                    ax[i,1].set_xticklabels(df.index)
                    ax[i,1].set_xlabel('Slab thickness / Å')
                    ax[i,1].set_ylabel('Time taken / s')
                    ax[i,1].plot(time, marker='x')
                    ax[i,1].legend(df.columns, title='Vacuum / Å')
                    ax[i,1].set_title('{}'.format(indices[i]))
                
                else: 
                    ax[i].set_prop_cycle(custom_cycler)
                    ax[i].set_xticks(list(range(len(df.index))))
                    ax[i].set_xticklabels(df.index)
                    ax[i].set_xlabel('Slab thickness / Å')
                    ax[i].set_ylabel('Surface energy / J m$^{-2}$')
                    ax[i].plot(val, marker='x')
                    ax[i].legend(df.columns, title='Vacuum / Å')
                    ax[i].set_title('{}'.format(indices[i]))
            
            plt.tight_layout()

    fig.savefig(plt_fname, bbox_inches='tight', facecolor='w')


def plot_enatom(df, time_taken=True, colors=None, dpi=300, width=6, height=5, 
heatmap=False, cmap='Wistia', plt_fname='energy_per_atom.png'):
    """
    Plots the energy per atom for all terminations. Based on surfaxe.convergence 
    parse_fols.

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols`, or any other 
            Dataframe with headings 'slab_thickness, 'vac_thickness', 
            'slab_per_atom', 'time_taken', 'index'. 
        time_taken (bool): Show the time taken for calculation to finish on the
            figure. Defaults to True.
        colors (`list`, optional): A list of colours for plots of different  
            vacuum thicknesses. Defaults to ``None``, which defaults to 
            surfaxe base style. 
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        width (`float`, optional): Width of figure in inches. Defaults to ``6``. 
        height (`float`, optional): Height of figure in inches. Defaults to 
            ``5``. 
        heatmap (`bool`, optional): If True plots a heatmap of surface energies.
            Defaults to False.
        cmap (`str`, optional): Matplotlib colourmap. Defaults to 'Wistia'
        plt_fname (`str`, optional): The name of the plot. Defaults to 
            ``energy_per_atom.png``. If name with no format suffix is supplied,  
            the format defaults to png.

    Returns:
        None, saves energy_per_atom.png
    """
    if not plt_fname.endswith('.png'):
        plt_fname += '.png'

    indices, vals, times, dfs, dfs_times = ([] for i in range(5))

    # Group the values by termination slab index, create df for time and
    # energy values. Converts the energy and time values to np arrays for 
    # plotting
    for group in df.groupby('slab_index'):
        df2 = group[1].pivot('slab_thickness', 'vac_thickness', 'slab_per_atom')
        df3 = group[1].pivot('slab_thickness', 'vac_thickness', 'time_taken')
        indices.append(group[0])
        vals.append(df2.to_numpy())
        times.append(df3.to_numpy())
        dfs.append(df2)
        dfs_times.append(df3)

    # Plots the heatmap
    if heatmap:
        if len(indices) == 1:
            fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))

            # Iterate through the values for plotting, create each plot on a 
            # separate ax, add the colourbar to each ax
            for index, val, time, df in zip(indices, vals, times, dfs):
                ax.set_yticks(list(range(len(df.index))))
                ax.set_yticklabels(df.index)
                ax.set_ylabel('Slab thickness / Å')
                ax.set_xticks(list(range(len(df.columns))))
                ax.set_xticklabels(df.columns)
                ax.set_xlabel('Vacuum thickness / Å')
                im = ax.imshow(val, cmap=cmap)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.2)
                cbar = plt.colorbar(im, cax=cax)
                cbar.set_label('Energy per atom / eV')
                ax.invert_yaxis()

            # Add the surface energy value labels to the plot, the for loop have 
            # to be in this order becuase it breaks the text if i and val are in 
            # the first loop, also j and k can't be zipped
            for df in dfs:
                for j in range(len(df.index)):
                    for k in range(len(df.columns)):
                        for val in vals:
                            ax.text(k, j, f"{val[j, k]: .3f}", ha="center",
                            va="bottom", color="black")

            # Add the time taken labels to the plot, same loop comment as above
            if time_taken:
                for df in dfs:
                    for j in range(len(df.index)):
                        for k in range(len(df.columns)):
                            for time in times:
                                ax.text(k, j, (f"{time[j, k]: .0f}"+' s'),
                                ha="center", va="top", color="black")

        # Plotting for multiple indices
        else:
            fig, ax = plt.subplots(ncols=len(indices), dpi=dpi, 
            figsize=(width, height))

            # Iterate through the values for plotting, create each plot on a 
            # separate ax, add the colourbar to each ax
            for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):
                ax[i].set_title('{}'.format(index))
                ax[i].set_yticks(list(range(len(df.index))))
                ax[i].set_yticklabels(df.index)
                ax[i].set_ylabel('Slab thickness / Å')
                ax[i].set_xticks(list(range(len(df.columns))))
                ax[i].set_xticklabels(df.columns)
                ax[i].set_xlabel('Vacuum thickness / Å')
                im = ax[i].imshow(val, cmap=cmap)
                divider = make_axes_locatable(ax[i])
                cax = divider.append_axes("right", size="5%", pad=0.2)
                cbar = plt.colorbar(im, cax=cax)
                cbar.set_label('Energy per atom / eV')
                ax[i].invert_yaxis()
            fig.tight_layout()

            # Add the surface energy value labels to the plot, the for loop have 
            # to be in this order becuase it breaks the text if i and val are in 
            # the first loop, also j and k can't be zipped
            for df in dfs:
                for j in range(len(df.index)):
                    for k in range(len(df.columns)):
                        for i, val in enumerate(vals):
                            ax[i].text(k, j, f"{val[j, k]: .3f}", ha="center",
                            va="bottom", color="black")

            # Add the time taken labels to the plot, same loop comment as above
            if time_taken:
                for df in dfs:
                    for j in range(len(df.index)):
                        for k in range(len(df.columns)):
                            for i, time in enumerate(times):
                                ax[i].text(k, j, (f"{time[j, k]: .0f}"+' s'),
                                ha="center", va="top", color="black")
    
    # Line plots
    else: 
        # Get the number of subplots 
        nrows = len(indices)
        ncols = 1
        if time_taken: 
            ncols = 2
        
        # Make the colour cycler with custom colours
        if colors is None:
            custom_cycler = (cycler(color=['#FFADB6','#FF596A','#CC4755',
                '#802D35','#E6505F','#40161A','#CC858C']))
        else: 
            custom_cycler = (cycler(color=colors))
        
        # Plot only the energy per atom for the only termination present 
        if nrows==1 and ncols==1: 
            fig, ax = plt.subplots(1,1, dpi=dpi, figsize=(width, height))
            for index, val, time, df in zip(indices, vals, times, dfs):
                ax.set_prop_cycle(custom_cycler)
                ax.set_title(index)
                ax.set_xticks(list(range(len(df.index))))
                ax.set_xticklabels(df.index)
                ax.set_xlabel('Slab thickness / Å')
                ax.set_ylabel('Energy per atom / eV')
                ax.plot(val, marker='x')
                ax.legend(df.columns, title='Vacuum / Å')

        # Plot energy per atom and time taken for the only termination present
        elif nrows==1 and ncols==2: 
            fig, ax = plt.subplots(1, 2, dpi=dpi, figsize=(width, height))
            for index, val, time, df in zip(indices, vals, times, dfs):
                ax[0].set_prop_cycle(custom_cycler)
                ax[0].set_title(index)
                ax[0].set_xticks(list(range(len(df.index))))
                ax[0].set_xticklabels(df.index)
                ax[0].set_xlabel('Slab thickness / Å')
                ax[0].set_ylabel('Energy per atom / eV')
                ax[0].plot(val, marker='x')
                ax[0].legend(df.columns, title='Vacuum / Å')
                
                ax[1].set_prop_cycle(custom_cycler)
                ax[1].set_title(index)
                ax[1].set_xticks(list(range(len(df.index))))
                ax[1].set_xticklabels(df.index)
                ax[1].set_xlabel('Slab thickness / Å')
                ax[1].set_ylabel('Time taken / s')
                ax[1].plot(time, marker='x')
                ax[1].legend(df.columns, title='Vacuum / Å')
                plt.tight_layout()

        # Plot energy per atom and time taken for different terminations 
        else: 
            fig, ax = plt.subplots(nrows=nrows,ncols=ncols, dpi=dpi, 
            figsize=(width, height))
            for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):
                # Separate plotting of times and energies, energies in the first
                # column, times in the second. Annoyingly matplotlib doesn't 
                # like ax[i,0] if there isn't a second axis so need to separate 
                # it for time_taken True and False
            
                if time_taken: 
                    ax[i,0].set_prop_cycle(custom_cycler)
                    ax[i,0].set_xticks(list(range(len(df.index))))
                    ax[i,0].set_xticklabels(df.index)
                    ax[i,0].set_xlabel('Slab thickness / Å')
                    ax[i,0].set_ylabel('Energy per atom / eV')
                    ax[i,0].plot(val, marker='x')
                    ax[i,0].legend(df.columns, title='Vacuum / Å')
                    ax[i,0].set_title('{}'.format(indices[i]))
             
                    ax[i,1].set_prop_cycle(custom_cycler)
                    ax[i,1].set_xticks(list(range(len(df.index))))
                    ax[i,1].set_xticklabels(df.index)
                    ax[i,1].set_xlabel('Slab thickness/ Å')
                    ax[i,1].set_ylabel('Time taken / s')
                    ax[i,1].plot(time, marker='x')
                    ax[i,1].legend(df.columns, title='Vacuum / Å')
                    ax[i,1].set_title('{}'.format(indices[i]))
                
                else: 
                    ax[i].set_prop_cycle(custom_cycler)
                    ax[i].set_xticks(list(range(len(df.index))))
                    ax[i].set_xticklabels(df.index)
                    ax[i].set_xlabel('Slab thickness / Å')
                    ax[i].set_ylabel('Energy per atom / eV')
                    ax[i].plot(val, marker='x')
                    ax[i].legend(df.columns, title='Vacuum / Å')
                    ax[i].set_title('{}'.format(indices[i]))
            
            plt.tight_layout()

    fig.savefig(plt_fname, bbox_inches='tight', facecolor='w')
