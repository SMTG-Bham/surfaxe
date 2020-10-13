# pymatgen
from pymatgen.io.vasp.sets import DictSet
from pymatgen import Structure

# Misc 
import pandas as pd 
import numpy as np 
import os
import warnings 

# Monkeypatching straight from Stackoverflow
def custom_formatwarning(message, category, filename, lineno, line=''):
    # Ignore everything except the message
    return 'UserWarning: ' + str(message) + '\n'

# Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
mpl.rcParams['figure.figsize'] = (10.0,8.0)
mpl.rcParams.update({'font.size': 14})

def save_slabs(list_of_slabs, **kwargs): 
    struc = Structure.from_file(filename=kwargs.get('structure'))
    bulk_name = struc.formula.replace(" ", "")

    if kwargs.get('make_fols') or kwargs.get('make_input_files'): 
        for slab in list_of_slabs:
            os.makedirs(os.path.join(os.getcwd(), r'{}/{}_{}_{}'.format(slab['hkl'],
            slab['slab_t'], slab['vac_t'], slab['s_index'])), exist_ok=True)

            # Makes all input files (KPOINTS, POTCAR, INCAR) based on the config
            # dictionary
            if kwargs.get('make_input_files'):
                vis = DictSet(structure=slab['slab'], **kwargs)
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

def plot_bond_analysis(df=None, filename=None, **kwargs): 
    """
    Plots the bond distance with respect to fractional coordinate. Used in 
    conjunction with surfaxe.analysis.bond_analysis.   

    Args:
        df (pandas DataFrame): DataFrame from surfaxe.analysis.bond_analysis
        filename (str): path to csv file with data from bond_analysis

    Returns:
        Plot
    """ 
    
    if filename: 
        df = pd.read_csv(filename)
    elif df: 
        df = df
    else: 
        warnings.formatwarning = custom_formatwarning
        warnings.warn('Data not supplied')

    colors = plt.rcParams["axes.prop_cycle"]()
    fig, axs = plt.subplots(nrows=len(kwargs.get('atoms')))

    # Iterate though the atom combinations, plot the graphs
    i=0
    for atom1, atom2 in kwargs.get('atoms'):
        c = next(colors)["color"]
        x = df['{}_c_coord'.format(atom1)]
        y = df['{}-{}_bond_distance'.format(atom1,atom2)]
        axs[i].scatter(x, y, marker='x', color=c)
        axs[i].set_ylabel("Bond distance / Å ")
        axs[i].legend(['{}-{} bond'.format(atom1, atom2)])
        i+=1
    
    plt.xlabel("Fractional coordinates")
    plt.savefig(**kwargs)
    

def plot_electrostatic_potential(df=None, filename=None, **kwargs): 
    """
    Plots the electrostatic potential along one direction. 

    Args: 
        df: pandas DataFrame from surfaxe.analysis.electrostatic_potential 
        filename (str): the filename of csv with potential

    """
    if df is not None: 
        df = df
    elif filename is not None: 
        df = pd.read_csv(filename)
    else: 
        warnings.formatwarning = custom_formatwarning
        warnings.warn('Data not supplied')

    # Plot both planar and macroscopic, save figure
    fig, ax = plt.subplots()
    ax.plot(df['planar'], label='planar')
    ax.plot(df['macroscopic'], label='macroscopic')
    ax.legend()
    plt.ylabel('Potential / eV')
    plt.savefig(**kwargs)

def plot_surfen(df, time_taken=True, cmap='Wistia', fmt='png', dpi=300, **kwargs):
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
    hkl = kwargs.get('hkl')
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
    dpi=dpi, bbox_inches='tight', **kwargs)


def plot_enatom(df, time_taken=True, cmap='Wistia', fmt='png', dpi=300, **kwargs):
    """
    Plots the energy per atom for all terminations. Based on surfaxe.convergence 
    parse_fols.

    Args:
        df (pandas DataFrame): DataFrame from `parse_fols` 
        time_taken (bool): whether it shows the time taken for calculation to
        finish on the graph; default=True
        cmap (str): Matplotlib colourmap used; defaut='Wistia'
        fmt (str): format for the output file; default='png'
        dpi (int): dots per inch; default=300

    Returns:
        hkl_energy_per_atom.png
    """
    hkl = kwargs.get('hkl')
    hkl_string = ''.join(map(str, hkl))

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


    plt.savefig('{}_energy_per_atom.{}'.format(hkl_string,fmt),
    dpi=dpi, bbox_inches='tight', **kwargs)
