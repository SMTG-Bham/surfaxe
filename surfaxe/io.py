# pymatgen
from pymatgen.io.vasp.sets import DictSet
from pymatgen import Structure

# Misc 
import pandas as pd 
import numpy as np 
import os
import warnings 
warnings.filterwarnings('once')
import json

# Monkeypatching straight from Stackoverflow
def _custom_formatwarning(message, category, filename, lineno, line=''):
    # Ignore everything except the message
    return 'UserWarning: ' + str(message) + '\n'

# Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
              ``_config_directory`` folder. If the filename does not exist,  
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
                    would be: ``MgO/POSCAR_001_20_30_1``
 
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
            os.makedirs(os.path.join(os.getcwd(), r'{}/{}_{}_{}'.format(slab['hkl'],
            slab['slab_t'], slab['vac_t'], slab['s_index'])), exist_ok=True)

            # Makes all input files (KPOINTS, POTCAR, INCAR) based on the config
            # dictionary
            if make_input_files:
                cd = load_config_dict(config_dict)
                vis = DictSet(slab['slab'], cd, **save_slabs_kwargs)
                vis.write_input(
                    r'{}/{}_{}_{}'.format(slab['hkl'], slab['slab_t'], 
                        slab['vac_t'],slab['s_index'])
                    )

            # Just makes the folders with POSCARs
            else:
                slab['slab'].to(fmt=fmt,
                filename=r'{}/{}_{}_{}/{}'.format(slab['hkl'],
                slab['slab_t'], slab['vac_t'], slab['s_index'], name))

    # Makes POSCAR_hkl_slab_vac_index files in the bulk_name folder
    else:
        os.makedirs(os.path.join(os.getcwd(), r'{}'.format(bulk_name)),
        exist_ok=True)
        for slab in list_of_slabs:
            slab['slab'].to(fmt=fmt,
            filename=r'{}/{}_{}_{}_{}_{}.vasp'.format(bulk_name, name, 
            slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index']))

def plot_bond_analysis(bonds, df=None, filename=None, 
plt_fname='bond_analysis.png', dpi=300): 
    """
    Plots the bond distance with respect to fractional coordinate. Used in 
    conjunction with surfaxe.analysis.bond_analysis.   

    Args:
        bonds (`list` of `tuples): List of bonds to compare
            e.g. [('Y', 'O'), ('Ti', 'S')]; order of bond pairs must be the same  
            as in the Dataframe or provided file.
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
    
    if filename is not None: 
        df = pd.read_csv(filename)
    elif df is not None: 
        df = df
    else: 
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Data not supplied')

    colors = plt.rcParams["axes.prop_cycle"]()
    fig, axs = plt.subplots(nrows=len(bonds))

    if len(bonds)==1: 
        for atom1, atom2 in bonds: 
            x = df['{}_c_coord'.format(atom1)]
            y = df['{}-{}_bond_distance'.format(atom1,atom2)]
            axs.scatter(x, y, marker='x')
            axs.set_ylabel("Bond distance / Å")
            axs.legend(['{}-{} bond'.format(atom1, atom2)])
    else: 
        # Iterate though the atom combinations, plot the graphs
        for i, (atom1, atom2) in enumerate(bonds):
            c = next(colors)["color"]
            x = df['{}_c_coord'.format(atom1)]
            y = df['{}-{}_bond_distance'.format(atom1,atom2)]
            axs[i].scatter(x, y, marker='x', color=c)
            axs[i].set_ylabel("Bond distance / Å")
            axs[i].legend(['{}-{} bond'.format(atom1, atom2)])
    
    plt.xlabel("Fractional coordinate in c")
    plt.savefig(plt_fname, dpi=dpi)
    

def plot_electrostatic_potential(df=None, filename=None, dpi=300, 
plt_fname='potential.png'): 
    """
    Plots the planar and macroscopic electrostatic potential along one 
    direction. Can take either a DataFrame or a potential.csv file as input. 

    Args: 
        df (`pandas DataFrame`, optional): pandas DataFrame from 
            surfaxe.analysis.electrostatic_potential. Defaults to ``None``. 
        filename (`str`, optional): The filename of csv file with potential 
            data. Defaults to ``None``. 
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'potential.png'``.
    
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
    plt.savefig(plt_fname, dpi=dpi)

def plot_surfen(df, hkl, time_taken=True, cmap='Wistia', fmt='png', dpi=300, 
heatmap=False):
    """
    Plots the surface energy for all terminations. Based on surfaxe.convergence 
    parse_fols. 

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols`, or any other 
            Dataframe with headings 'slab_thickness, 'vac_thickness', 
            'surface_energy', 'time_taken', 'index'. 
        time_taken (bool): Show the time taken for calculation to finish on the
            figure. Defaults to True.
        cmap (`str`, optional): Matplotlib colourmap. Defaults to 'Wistia'
        fmt (`str`, optional): Format for the output file. Defaults to 'png'
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        heatmap (`bool`, optional): If True plots a heatmap of surface energies.
            Defaults to False.

    Returns:
        hkl_surface_energy.png
    """
    hkl_string = ''.join(map(str, hkl))

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
            fig, ax = plt.subplots(1,1)

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
                cbar = fig.colorbar(im, cax=cax, orientation='vertical')
                cbar.set_title('Surface energy / Jm$^{-2}$')
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
            mpl.rcParams['figure.figsize'] = (10.0,10.0)
            fig, ax = plt.subplots(ncols=len(indices))

            # Iterate through the values for plotting, create each plot on a 
            # separate ax, add the colourbar to each ax
            for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):
                ax[i].set_title('{}'.format(index))
                ax[i].set_yticks(list(range(len(df.index))))
                ax[i].set_yticklabels(df.index)
                ax[i].set_ylabel('Slab thickness')
                ax[i].set_xticks(list(range(len(df.columns))))
                ax[i].set_xticklabels(df.columns)
                ax[i].set_xlabel('Vacuum thickness')
                im = ax[i].imshow(val, cmap=cmap)
                divider = make_axes_locatable(ax[i])
                cax = divider.append_axes("right", size="5%", pad=0.2)
                cbar = plt.colorbar(im, cax=cax)
                cbar.set_title('Surface energy / Jm$^{-2}$')
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
        # Get the number of subplots 
        nrows = len(indices)
        ncols = 1
        if time_taken: 
            ncols = 2
        
        # Plot only the surface energy for the only termination present 
        if nrows==1 and ncols==1: 
            fig, ax = plt.subplots(1,1)
            dfs[0].plot(ax=ax, marker='x', xticks=dfs[0].index, 
            xlabel='Slab thickness / Å', ylabel='Surface energy / Jm$^{-2}$')
            ax.legend(title='Vacuum / Å')

        elif nrows==1 and ncols==2: 
            fig, ax = plt.subplots(1,2)
            dfs[0].plot(ax=ax[0], xticks=dfs[0].index, marker='x', 
            xlabel='Slab thickness / Å', ylabel='Surface energy / Jm$^{-2}$')
            ax[0].legend(title='Vacuum / Å')
            dfs_times[0].plot(ax=ax[1], xticks=dfs_times[0].index, marker='x', 
            xlabel='Slab thickness / Å', ylabel='Time taken / s')
            ax[1].legend(title='Vacuum / Å')
            plt.tight_layout()

        # Plot times taken and or different terminations 
        else: 
            fig, ax = plt.subplots(nrows=nrows,ncols=ncols) 

            # Separate plotting of times and energies, energies in the first
            # column, times in the second. Annoyingly matplotlib doesn't like 
            # ax[i,0] if there isn't a second axis so need to separate it for 
            # time_take True and False
            
            if time_taken: 
                for i, df in enumerate(dfs): 
                    df.plot(ax=ax[i,0], xticks=df.index, marker='x', 
                    xlabel='Slab thickness / Å', 
                    ylabel='Surface energy / Jm$^{-2}$')
                    ax[i,0].legend(title='Vacuum / Å')
                    ax[i,0].set_title('{}'.format(indices[i]))
            
                for i, df in enumerate(dfs_times): 
                    df.plot(ax=ax[i,1], xticks=df.index, marker='x', 
                    xlabel='Slab thickness / Å', ylabel='Time taken / s')
                    ax[i,1].legend(title='Vacuum / Å')
                    ax[i,1].set_title('{}'.format(indices[i]))
            else: 
                for i, df in enumerate(dfs): 
                    df.plot(ax=ax[i], xticks=df.index, marker='x', 
                    xlabel='Slab thickness / Å', 
                    ylabel='Surface energy / Jm$^{-2}$')
                    ax[i].legend(title='Vacuum / Å')
                    ax[i].set_title('{}'.format(indices[i]))
            
            plt.tight_layout()

    plt.savefig('{}_surface_energy.{}'.format(hkl_string, fmt),
    dpi=dpi, bbox_inches='tight')


def plot_enatom(df, hkl, time_taken=True, cmap='Wistia', fmt='png', dpi=300, 
heatmap=False):
    """
    Plots the energy per atom for all terminations. Based on surfaxe.convergence 
    parse_fols.

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols`, or any other 
            Dataframe with headings 'slab_thickness, 'vac_thickness', 
            'slab_per_atom', 'time_taken', 'index'. 
        time_taken (bool): Show the time taken for calculation to finish on the
            figure. Defaults to True.
        cmap (`str`, optional): Matplotlib colourmap. Defaults to 'Wistia'
        fmt (`str`, optional): Format for the output file. Defaults to 'png'
        dpi (`int`, optional): Dots per inch. Defaults to 300.
        heatmap (`bool`, optional): If True plots a heatmap of surface energies.
            Defaults to False.

    Returns:
        hkl_energy_per_atom.png
    """

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
            fig, ax = plt.subplots(1,1)

            # Iterate through the values for plotting, create each plot on a 
            # separate ax, add the colourbar to each ax
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
            fig, ax = plt.subplots(ncols=len(indices))

            # Iterate through the values for plotting, create each plot on a 
            # separate ax, add the colourbar to each ax
            for i, (index, val, time, df) in enumerate(zip(indices, vals, times, dfs)):
                ax[i].set_title('{}'.format(index))
                ax[i].set_yticks(list(range(len(df.index))))
                ax[i].set_yticklabels(df.index)
                ax[i].set_ylabel('Slab thickness')
                ax[i].set_xticks(list(range(len(df.columns))))
                ax[i].set_xticklabels(df.columns)
                ax[i].set_xlabel('Vacuum thickness')
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
        
        # Plot only the surface energy for the only termination present 
        if nrows==1 and ncols==1: 
            fig, ax = plt.subplots(1,1)
            dfs[0].plot(ax=ax, xticks=df.index, marker='x', 
            xlabel='Slab thickness / Å', ylabel='Energy per atom / eV')
            ax.legend(title='Vacuum / Å')

        elif nrows==1 and ncols==2: 
            fig, ax = plt.subplots(1,2)
            dfs[0].plot(ax=ax[0], xticks=dfs[0].index, marker='x', 
            xlabel='Slab thickness / Å', ylabel='Energy per atom / eV')
            ax[0].legend(title='Vacuum / Å')
            dfs_times[0].plot(ax=ax[1], xticks=dfs_times[0].index, marker='x', 
            xlabel='Slab thickness / Å', ylabel='Time taken / s')
            ax[1].legend(title='Vacuum / Å')
            plt.tight_layout()

        # Plot times taken and or different terminations 
        else: 
            fig, ax = plt.subplots(nrows=nrows,ncols=ncols)

            # Separate plotting of times and energies, energies in the first
            # column, times in the second. Annoyingly matplotlib doesn't like 
            # ax[i,0] if there isn't a second axis so need to separate it for 
            # time_take True and False
            
            if time_taken: 
                for i, df in enumerate(dfs): 
                    df.plot(ax=ax[i,0], xticks=df.index, marker='x', 
                    xlabel='Slab thickness / Å', ylabel='Energy per atom / eV')
                    ax[i,0].legend(title='Vacuum / Å')
                    ax[i,0].set_title('{}'.format(indices[i]))
            
                for i, df in enumerate(dfs_times): 
                    df.plot(ax=ax[i,1], xticks=df.index, marker='x', 
                    xlabel='Slab thickness / Å', ylabel='Time taken / s')
                    ax[i,1].legend(title='Vacuum / Å')
                    ax[i,1].set_title('{}'.format(indices[i]))
            else: 
                for i, df in enumerate(dfs): 
                    df.plot(ax=ax[i], xticks=df.index, marker='x', 
                    xlabel='Slab thickness / Å', ylabel='Energy per atom / eV')
                    ax[i].legend(title='Vacuum / Å')
                    ax[i].set_title('{}'.format(indices[i]))
            
            plt.tight_layout()

    hkl_string = ''.join(map(str, hkl))
    plt.savefig('{}_energy_per_atom.{}'.format(hkl_string,fmt),
    dpi=dpi, bbox_inches='tight')
