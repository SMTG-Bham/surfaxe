# pymatgen
from pymatgen.io.vasp.sets import DictSet
from pymatgen import Structure

# Misc 
import pandas as pd 
import numpy as np 
import os
import warnings 
import json

# Monkeypatching straight from Stackoverflow
def _custom_formatwarning(message, category, filename, lineno, line=''):
    # Ignore everything except the message
    return 'UserWarning: ' + str(message) + '\n'

# Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
mpl.rcParams['figure.figsize'] = (10.0,8.0)
mpl.rcParams.update({'font.size': 14})

from surfaxe import _config_directory

def load_config_dict(config_dict): 
    """
    Loads the config dictionary for writing VASP input files. 

    Args: 
        config_dict(``None``, `dict` or `str`): The config dict containing info  
            on INCAR, POTCAR and KPOINTS settings. Can be supplied as: 

            * ``dict``: All settings for the calculations provided as a 
              dictionary of dictionaries 

                    e.g. ``{'INCAR': {'ENCUT': 500, 'ISYM': 2, 'GGA': 'PE'}, 
                    'KPOINTS': {'reciprocal_density': 20}, 
                    'POTCAR': {'Sn': 'Sn_d', 'O': 'O'}}``

            * ``str``: Filename of the config dictionary in the 
              ``_config_dictionaries`` folder. If the filename does not exist,  
              the function defaults to the ``PBEsol_config.json`` file.            

            * ``None``: The default option, makes a PBEsol config dictionary for
              a single shot calculation from the ``PBEsol_config.json`` file.  
    Returns: 
        Dictionary
    
    """
    if type(config_dict) is dict: 
        cd = config_dict 
    elif type(config_dict) is str: 
        if os.path.isfile(os.path.join(_config_directory, config_dict)): 
            with open(os.path.join(_config_directory, config_dict), 'r') as f:
                cd = json.load(f)
        else: 
            with open(os.path.join(_config_directory, 'PBEsol_config.json'), 'r') as f:
                cd = json.load(f)
    else: 
        with open(os.path.join(_config_directory, 'PBEsol_config.json'), 'r') as f:
            cd = json.load(f)

    return cd 

def slabs_to_file(list_of_slabs, structure, make_fols, make_input_files, 
config_dict, potcar_functional, user_incar_settings, user_kpoints_settings, 
user_potcar_settings): 
    """
    Saves the slabs to file 

    Args: 
        list_of_slabs (list): a list of slab dictionaries made with either of
        surfaxe.generation get_slab functions 

    Returns: 
        POSCARs in folders 
    """
    struc = Structure.from_file(structure)
    bulk_name = struc.formula.replace(" ", "")

    if make_fols or make_input_files: 
        for slab in list_of_slabs:
            os.makedirs(os.path.join(os.getcwd(), r'{}/{}_{}_{}'.format(slab['hkl'],
            slab['slab_t'], slab['vac_t'], slab['s_index'])), exist_ok=True)

            # Makes all input files (KPOINTS, POTCAR, INCAR) based on the config
            # dictionary
            if make_input_files:
                cd = load_config_dict(config_dict)
                vis = DictSet(slab['slab'], cd,  
                user_incar_settings=user_incar_settings, 
                user_kpoints_settings=user_kpoints_settings, 
                user_potcar_settings=user_potcar_settings)
                vis.write_input(r'{}/{}_{}_{}'.format(slab['hkl'],
                                                      slab['slab_t'],
                                                      slab['vac_t'],
                                                      slab['s_index']))

            # Just makes the folders with POSCARs
            else:
                slab['slab'].to(fmt='poscar',
                filename=r'{}/{}_{}_{}/POSCAR'.format(slab['hkl'],
                slab['slab_t'], slab['vac_t'], slab['s_index']))

    # Makes POSCAR_hkl_slab_vac_index files in the bulk_name folder
    else:
        os.makedirs(os.path.join(os.getcwd(), r'{}'.format(bulk_name)),
        exist_ok=True)
        for slab in list_of_slabs:
            slab['slab'].to(fmt='poscar',
            filename=r'{}/POSCAR_{}_{}_{}_{}.vasp'.format(bulk_name,slab['hkl'],
            slab['slab_t'], slab['vac_t'], slab['s_index']))

def plot_bond_analysis(bonds, df=None, filename=None, plt_fname='bond_analysis.png', 
dpi=300): 
    """
    Plots the bond distance with respect to fractional coordinate. Used in 
    conjunction with surfaxe.analysis.bond_analysis.   

    Args:
        bonds (`list` of `tuples): List of bonds to compare
            e.g. [('Y', 'O'), ('Ti', 'S')]; order does not matter
        df (`pandas DataFrame`, optional): DataFrame from 
            surfaxe.analysis.bond_analysis. Defaults to ``None``.
        filename (`str`, optional): Path to csv file with data from 
            surfaxe.analysis.bond_analysis. Defaults to ``None``. 
            Either df or filename need to be supplied.
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'bond_analysis.png'``.
        dpi (`int`, optional): Dots per inch. Defaults to ``300``.  

    Returns:
        Plot
    """ 
    
    if filename: 
        df = pd.read_csv(filename)
    elif df: 
        df = df
    else: 
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Data not supplied')

    colors = plt.rcParams["axes.prop_cycle"]()
    fig, axs = plt.subplots(nrows=len(bonds))

    # Iterate though the atom combinations, plot the graphs
    i=0
    for atom1, atom2 in bonds:
        c = next(colors)["color"]
        x = df['{}_c_coord'.format(atom1)]
        y = df['{}-{}_bond_distance'.format(atom1,atom2)]
        axs[i].scatter(x, y, marker='x', color=c)
        axs[i].set_ylabel("Bond distance / Å ")
        axs[i].legend(['{}-{} bond'.format(atom1, atom2)])
        i+=1
    
    plt.xlabel("Fractional coordinates")
    plt.savefig(plt_fname, dpi)
    

def plot_electrostatic_potential(df=None, filename=None, dpi=300,
    plt_fname='electrostatic_potential.png'): 
    """
    Plots the electrostatic potential along one direction. 

    Args: 
        df: pandas DataFrame from surfaxe.analysis.electrostatic_potential 
        filename (str): the filename of csv with potential
    
    Returns: 
        potential.png
    """
    if df is not None: 
        df = df
    elif filename is not None: 
        df = pd.read_csv(filename)
    else: 
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Data not supplied')

    # Plot both planar and macroscopic, save figure
    fig, ax = plt.subplots()
    ax.plot(df['planar'], label='planar')
    ax.plot(df['macroscopic'], label='macroscopic')
    ax.legend()
    plt.ylabel('Potential / eV')
    plt.savefig(plt_fname, dpi)

def plot_surfen(df, hkl=None, time_taken=True, cmap='Wistia', fmt='png', dpi=300):
    """
    Plots the surface energy for all terminations. Based on surfaxe.convergence 
    parse_fols. 

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols`
        time_taken (bool): whether it shows the time taken for calculation to
            finish on the graph; default=True
        cmap (str): Matplotlib colourmap; defaut='Wistia'
        fmt (str): format for the output file; default='png'
        dpi (int): dots per inch; default=300

    Returns:
        hkl_surface_energy.png
    """
    hkl_string = ''.join(map(str, hkl))

    indices, vals, times, dfs = ([] for i in range(4))

    # Group the values by termination slab index, create df for time and
    # energy values. Converts the energy and time values to np arrays for plotting
    for group in df.groupby('slab_index'):
        df2 = group[1].pivot(index='slab_thickness',
                             columns='vac_thickness',
                             values='surface_energy')
        df3 = group[1].pivot(index='slab_thickness',
                             columns='vac_thickness',
                             values='time_taken')
        indices.append(group[0])
        vals.append(df2.to_numpy())
        times.append(df3.to_numpy())
        dfs.append(df2)

    # Plotting has to be separated into plotting for one and more than one
    # indices mpl won't let you index axes if there's only one set
    if len(indices) == 1:
        mpl.rcParams['figure.figsize'] = (6.0,6.0)
        fig, ax = plt.subplots(1,1)
        ax.set_title('{} surface energies'.format(hkl))

        for (index, val, time, df) in zip(indices, vals, times, dfs):
            ax.set_yticks(list(range(len(df.index))))
            ax.set_yticklabels(df.index)
            ax.set_ylabel('Slab thickness')
            ax.set_xticks(list(range(len(df.columns))))
            ax.set_xticklabels(df.columns)
            ax.set_xlabel('Vacuum thickness')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            im = ax.imshow(val, cmap=cmap)
            fig.colorbar(im, cax=cax, orientation='vertical')
            ax.invert_yaxis()

        # Add the surface energy value labels to the plot - the for loops are
        # needed in this order and they can't be zipped because this is the
        # only way they display the values correctly
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
        mpl.rcParams['figure.figsize'] = (10.0,10.0)
        fig, ax = plt.subplots(ncols=len(indices))

        # The extra args are there because the default leaves a massive gap
        # between the title and subplot titles, no amount of changing tight_layout
        # and subplots_adjust helped with the issue - this is a known mpl issue
        fig.suptitle('{} surface energies'.format(hkl), y=0.73, ha='center',
                     va='top', size=22)

        # Iterate through the values for plotting, create each plot on a separate 
        # ax, add the colourbar to each ax
        for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):
            ax[i].set_title('Slab index {}'.format(index))
            ax[i].set_yticks(list(range(len(df.index))))
            ax[i].set_yticklabels(df.index)
            ax[i].set_ylabel('Slab thickness')
            ax[i].set_xticks(list(range(len(df.columns))))
            ax[i].set_xticklabels(df.columns)
            ax[i].set_xlabel('Vacuum thickness')
            im = ax[i].imshow(val, cmap=cmap)
            divider = make_axes_locatable(ax[i])
            cax = divider.append_axes("right", size="5%", pad=0.2)
            plt.colorbar(im, cax=cax)
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


    plt.savefig('{}_surface_energy.{}'.format(hkl_string, fmt),
    dpi=dpi, bbox_inches='tight')


def plot_enatom(df, hkl=None, time_taken=True, cmap='Wistia', fmt='png', dpi=300):
    """
    Plots the energy per atom for all terminations. Based on surfaxe.convergence 
    parse_fols.

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols` 
        time_taken (`bool`): whether it shows the time taken for calculation to
            finish on the graph; default=True
        cmap (`str`): Matplotlib colourmap used; defaut='Wistia'
        fmt (`str`): format for the output file; default='png'
        dpi (`int`): dots per inch; default=300

    Returns:
        hkl_energy_per_atom.png
    """

    indices, vals, times, dfs = ([] for i in range(4))

    # Group the values by termination slab index, create df for time and
    # energy values. Converts the energy and time values to np arrays for plotting
    for group in df.groupby('slab_index'):
        df2 = group[1].pivot(index='slab_thickness',
                             columns='vac_thickness',
                             values='slab_per_atom')
        df3 = group[1].pivot(index='slab_thickness',
                             columns='vac_thickness',
                             values='time_taken')
        indices.append(group[0])
        vals.append(df2.to_numpy())
        times.append(df3.to_numpy())
        dfs.append(df2)

    if len(indices) == 1:
        mpl.rcParams['figure.figsize'] = (6.0,6.0)
        fig, ax = plt.subplots(1,1)
        ax.set_title('{} energies per atom'.format(hkl))

        # Iterate through the values for plotting, create each plot on a separate ax,
        # add the colourbar to each ax
        for index, val, time, df in zip(indices, vals, times, dfs):
            ax.set_yticks(list(range(len(df.index))))
            ax.set_yticklabels(df.index)
            ax.set_ylabel('Slab thickness')
            ax.set_xticks(list(range(len(df.columns))))
            ax.set_xticklabels(df.columns)
            ax.set_xlabel('Vacuum thickness')
            im = ax.imshow(val, cmap=cmap)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            plt.colorbar(im, cax=cax)
            ax.invert_yaxis()

        # Add the surface energy value labels to the plot, the for loop have to
        # be in this order becuase it breaks the text if i and val are in the
        # first loop, also j and k can't be zipped
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
        mpl.rcParams['figure.figsize'] = (10.0,10.0)
        fig, ax = plt.subplots(ncols=len(indices))

        # The extra args are there because the default leaves a massive gap
        # between the title and subplot titles, no amount of changing tight_layout
        # and subplots_adjust helped with the issue
        fig.suptitle('{} energies per atom'.format(hkl), y=0.73, ha='center',
                     va='top', size=22)

        # Iterate through the values for plotting, create each plot on a separate ax,
        # add the colourbar to each ax
        for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):
            ax[i].set_title('Slab index {}'.format(index))
            ax[i].set_yticks(list(range(len(df.index))))
            ax[i].set_yticklabels(df.index)
            ax[i].set_ylabel('Slab thickness')
            ax[i].set_xticks(list(range(len(df.columns))))
            ax[i].set_xticklabels(df.columns)
            ax[i].set_xlabel('Vacuum thickness')
            im = ax[i].imshow(val, cmap=cmap)
            divider = make_axes_locatable(ax[i])
            cax = divider.append_axes("right", size="5%", pad=0.2)
            plt.colorbar(im, cax=cax)
            ax[i].invert_yaxis()
        fig.tight_layout()

        # Add the surface energy value labels to the plot, the for loop have to
        # be in this order becuase it breaks the text if i and val are in the
        # first loop, also j and k can't be zipped
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

    hkl_string = ''.join(map(str, hkl))
    plt.savefig('{}_energy_per_atom.{}'.format(hkl_string,fmt),
    dpi=dpi, bbox_inches='tight')
