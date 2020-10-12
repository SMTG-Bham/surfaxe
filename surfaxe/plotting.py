# pymatgen
from pymatgen.io.vasp.sets import DictSet
from pymatgen import Structure

# Misc 
import pandas as pd 
import numpy as np 
import os
import warnings 

# Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
mpl.rcParams['figure.figsize'] = (10.0,8.0)
mpl.rcParams.update({'font.size': 14})

# surfaxe 
from surfaxe.generation import get_all_slabs, get_one_hkl_slabs
from surfaxe.analysis import electrostatic_potential, bond_analysis, simple_nn,\
     complex_nn, cart_displacements
from surfaxe.convergence import parse_fols

# in theory should make this easier to use with python?
# in practice who knows 

def save_csv(func, *args, **kwargs): 
    df = func(*args, **kwargs)
    if 'csv_fname' in kwargs: 
        csv_fname = kwargs.get('csv_fname')
    else: 
        csv_fname = 'data.csv'
    df.to_csv(csv_fname, header=True, index=False, **kwargs)

def save_txt(func, *args, **kwargs): 
    df = func(*args, **kwargs)
    if 'txt_fname' in kwargs: 
        txt_fname = kwargs.get('txt_fname')
    else: 
        txt_fname = 'data.txt'
    df.to_csv(txt_fname, header=True, index=False, sep='\t', mode='w')

def save_slabs(func, *args, **kwargs): 
    big_list = func(*args, **kwargs)
    struc = Structure.from_file(filename=kwargs.get('structure'))
    bulk_name = struc.formula.replace(" ", "")

    if kwargs.get('make_fols') or kwargs.get('make_input_files'): 
        for slab in big_list:
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
        for slab in func(*args, **kwargs):
            slab['slab'].to(fmt='poscar',
            filename=r'{}/POSCAR_{}_{}_{}_{}.vasp'.format(bulk_name,slab['hkl'],
            slab['slab_t'], slab['vac_t'], slab['s_index']))

def plot_bond_analysis(bond_analysis, plt_fname='bond_analysis.png', dpi=300, **kwargs): 
    """
    Plots the bond distance with respect to fractional coordinate graph from the
    csv file generated with `bond_analysis`. The required argument is `atoms`. 

    Args:
        atoms (list of tuples) in the same order as in bond_analysis
        plt_fname (str): filename of the plot
        dpi (int): dots per inch; default=300

    Returns:
        Plot
    """ 
    # Check if a bond_analysis file exists, otherwise get the dataframe directly
    # from the function.
    if os.path.isfile('./bond_analysis.csv'): 
        df = pd.read_csv('./bond_analysis.csv')

    else: 
        df = bond_analysis(**kwargs)

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
    plt.savefig(plt_fname, dpi=dpi, **kwargs)
    

def plot_electrostatic_potential(*args, plt_fname='potential.png', dpi=300, 
**kwargs): 
    if os.path.isfile('./potential.csv'): 
        df = pd.read_csv('./potential.csv')

    else: 
        df = electrostatic_potential(**kwargs)
    
    # Plot both planar and macroscopic, save figure
    fig, ax = plt.subplots()
    ax.plot(df['planar'], label='planar')
    ax.plot(df['macroscopic'], label='macroscopic')
    ax.legend()
    plt.ylabel('Potential / eV')
    plt.savefig(plt_fname, dpi=dpi, **kwargs)

def plot_surfen(parse_fols, hkl=None, time_taken=True, cmap='Wistia', fmt='png', dpi=300, **kwargs):
    """
    Reads from the csv file created by `parse_fols` to plot the surface energy
    for all possible terminations

    Args:
        hkl (tuple): Miller index
        time_taken (bool): whether it shows the time taken for calculation to
        finish on the graph; default=True
        cmap (str): Matplotlib colourmap; defaut='Wistia'
        fmt (str): format for the output file; default='png'
        dpi (int): dots per inch; default=300

    Returns:
        hkl_surface_energy.png
    """
    # Check all neccessary input parameters are present 
    if not hkl: 
        raise ValueError('The required argument hkl was not supplied')

    hkl_string = ''.join(map(str, hkl))

    if os.path.isfile('./{}_data.csv'.format(hkl_string)): 
        df = pd.read_csv('./{}_data.csv'.format(hkl_string))

    else: 
        df = parse_fols(**kwargs)

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


    plt.savefig('{}_surface_energy.{}'.format(''.join(map(str, hkl)), fmt),
    dpi=dpi, bbox_inches='tight', **kwargs)


def plot_enatom(hkl=None, time_taken=True, cmap='Wistia', fmt='png', dpi=300, **kwargs):
    """
    Reads from the csv file created by `parse_fols` to plot the energy per atom
    for all possible terminations

    Args:
        hkl (tuple): Miller index
        time_taken (bool): whether it shows the time taken for calculation to
        finish on the graph; default=True
        cmap (str): Matplotlib colourmap used; defaut='Wistia'
        fmt (str): format for the output file; default='png'
        dpi (int): dots per inch; default=300

    Returns:
        hkl_energy_per_atom.png
    """
    # Check all neccessary input parameters are present 
    if not hkl: 
        raise ValueError('The required argument hkl was not supplied') 

    hkl_string = ''.join(map(str, hkl))


    df = pd.read_csv('{}_data.csv'.format(hkl_string))

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
