
from pymatgen import Structure, Element, Specie
from pymatgen.core.structure import SiteCollection
from pymatgen.core.lattice import Lattice
from pymatgen.analysis.local_env import BrunnerNN_real, BrunnerNN_reciprocal,\
BrunnerNN_relative, CovalentBondNN, Critic2NN, CrystalNN, CutOffDictNN, EconNN,\
JmolNN, MinimumDistanceNN, MinimumOKeeffeNN, MinimumVIRENN, NearNeighbors, VoronoiNN
from pymatgen.io.vasp.outputs import Locpot, VolumetricData
import os
import math
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('once')

## Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
mpl.rcParams['figure.figsize'] = (10.0,8.0)
mpl.rcParams.update({'font.size': 14})


def cart_displacements(start, end, elements, max_disp=0.1,
                       csv_fname='cart_displacements.txt'):
    """
    Produces a text file with all the magnitude of displacements of atoms
    in Cartesian space

    Args:
        start (str): filename of initial structure file in any format supported
        by pymatgen.
        end (str): filename of final structure file in any format supported
        by pymatgen.
        elements (list): list of elements in the structure
        e.g. ['Y', 'Ti', 'O', 'S'] in any order
        max_disp (float): maximum displacement shown; default 0.1 Å
        csv_fname (str): filename of the file produced

    Returns:
        Text file with displacements

    """
    # Instantiate the structures from files
    start_struc = Structure.from_file(start)
    end_struc = Structure.from_file(end)

    # Add the site labels to the structure
    el_dict = {i : 1 for i in elements}
    site_labels = []

    for site in start_struc:
        symbol = site.specie.symbol
        site_labels.append((symbol,el_dict[symbol]))
        el_dict[symbol] +=1
    start_struc.add_site_property('', site_labels)

    # Convert to cartesian coordinates
    start_struc = start_struc.cart_coords
    end_struc = end_struc.cart_coords

    # Calculate the displacements
    disp_list = []
    for n, (start_coord, end_coord) in enumerate(zip(start_struc, end_struc)):
        xdisp = math.pow(start_coord[0] - end_coord[0], 2)
        ydisp = math.pow(start_coord[1] - end_coord[1], 2)
        zdisp = math.pow(start_coord[2] - end_coord[2], 2)
        d = math.sqrt(xdisp + ydisp + zdisp)
        label = site_labels[n]
        if d >= max_disp:
            disp_list.append({'site': n+1,
                             'atom': label,
                             'displacement': f"{d: .3f}"})
    # Save as txt file
    df = pd.DataFrame(disp_list)
    df.to_csv(csv_fname, header=True, index=False, sep='\t', mode='w')

def bond_analysis(structure, atoms, nn_method=CrystalNN(), ox_states=None,
                  return_df=False, **kwargs):
    """
    Parses the structure looking for bonds between atoms. Check the validity of
    nearest neighbour method on the bulk structure before using it on slabs.

    Args:
        structure (str): filename of structure, takes all pymatgen-supported formats.
        atoms (list of tuples): list of bonds to compare
        e.g. [('Y', 'O'), ('Ti', 'S')]; order does not matter
        nn_method (class): pymatgen.analysis.local_env nearest neighbour method;
        default=CrystalNN()
        ox_states (list or dict): add oxidation states either by sites
        i.e. [3, 2, 2, 1, -2, -2, -2, -2] or by element i.e. {'Fe': 3, 'O':-2};
        default=None which adds oxidation states by guess
        return_df (bool): returns the DataFrame; default=False

    Returns:
        csv file
    """

    struc = Structure.from_file(structure)

    # Adds oxidation states by guess by default or if the provided oxidation states are
    # antyhing but a list or a dict; max_sites speeds up the by_guess method
    if type(ox_states) is dict:
        struc.add_oxidation_state_by_element(ox_states)
    elif type(ox_states) is list:
        struc.add_oxidation_state_by_site(ox_states)
    else:
        struc.add_oxidation_state_by_guess(max_sites=-1)

    # Iterates through the structure, looking for pairs of bonded atoms. If the
    # sites match, the bond distance is calculated and passed to a dataframe
    bonds_info = []
    for n, pos in enumerate(struc):
        for atom1, atom2 in atoms:
            if pos.specie.symbol is atom1:
                nearest_neighbours = nn_method.get_nn_info(struc, n)
                matched_sites = []
                for d in nearest_neighbours:
                    if d.get('site').specie.symbol is atom2:
                        matched_sites.append(d)
                bond_distances = [struc.get_distance(n,x['site_index']) for x in matched_sites]
                bonds_info.append({'{}_index'.format(atom1): n+1,
                                   '{}_c_coord'.format(atom1): pos.c,
                                   '{}-{}_bond_distance'.format(atom1,atom2): np.mean(bond_distances)})

    df = pd.DataFrame(bonds_info)
    df.to_csv('bond_analysis_data.csv', index=False)

    if return_df is True:
        return df

def plot_bond_analysis(atoms, plt_fname='bond_analysis.png', dpi=300, **kwargs):
    """
    Plots the bond distance with respect to fractional coordinate graph from the
    csv file generated with `bond_analysis`

    Args:
        atoms (list of tuples) in the same order as in bond_analysis
        plt_fname (str): filename of the plot
        dpi (int): dots per inch; default=300

    Returns:
        Plot
    """
    # Read in csv, define colour wheel
    df = pd.read_csv('bond_analysis_data.csv')
    colors = plt.rcParams["axes.prop_cycle"]()

    fig, axs = plt.subplots(nrows=len(atoms))

    # Iterate though the atom combinations, plot the graphs
    i=0
    for atom1, atom2 in atoms:
        c = next(colors)["color"]
        x = df['{}_c_coord'.format(atom1)]
        y = df['{}-{}_bond_distance'.format(atom1,atom2)]
        axs[i].scatter(x, y, marker='x', color=c)
        axs[i].set_ylabel("Bond distance / Å ")
        axs[i].legend(['{}-{} bond'.format(atom1, atom2)])
        i+=1
    plt.xlabel("Fractional coordinates")
    plt.savefig(plt_fname, dpi=dpi)

def electrostatic_potential(lattice_vector, filename='./LOCPOT', axis=2,
                            make_csv=True, csv_fname='potential.csv',
                            plt_fname='potential.png', dpi=300, **kwargs):
    """
    Reads LOCPOT to get the planar and macroscopic potential in specified direction

    Args:
        lattice_vector (float): the periodicity of the slab
        filename (str): path to your locpot file, default='./LOCPOT'
        axis (int): direction in which the potential is investigated; a=0, b=1,
        c=2; default=2
        make_csv (bool): makes a csv file with planar and macroscopic potential,
        default=True
        csv_fname (str): filename of the csv file, default='potential.csv'
        plt_fname (str): filename of the plot of potentials, controls the format,
        default='potential.png'
        dpi (int): dots per inch; default=300

    Returns:
        csv file and plot of planar and macroscopic potential
    """
    # Read potential and structure data
    lpt = Locpot.from_file(filename)
    struc = Structure.from_file(filename)

    # Planar potential
    planar = lpt.get_average_along_axis(axis)

    # Divide lattice parameter by no. of grid points in the direction
    resolution = struc.lattice.abc[axis]/lpt.dim[axis]

    # Get number of points over which the rolling average is evaluated
    points = int(lattice_vector/resolution)

    # Need extra points at the start and end of planar potential to evaluate the
    # macroscopic potential this makes use of the PBC where the end of one unit
    # cell coincides with start of the next one
    add_to_start = planar[(len(planar) - points): ]
    add_to_end = planar[0:points]
    pfm_data = np.concatenate((add_to_start,planar,add_to_end))
    pfm = pd.DataFrame(data=pfm_data, columns=['y'])

    # Macroscopic potential
    m_data = pfm.y.rolling(window=points, center=True).mean()
    macroscopic = m_data.iloc[points:(len(planar)+points)]
    macroscopic.reset_index(drop=True,inplace=True)

    # Make csv
    if make_csv is True:
        data = pd.DataFrame(data=planar, columns=['planar'])
        data['macroscopic'] = macroscopic
        data.to_csv(csv_fname, header = True, index = False)

    # Plot both planar and macroscopic, save figure
    fig,ax = plt.subplots()
    ax.plot(planar, label='planar')
    ax.plot(macroscopic, label='macroscopic')
    ax.legend()
    plt.ylabel('Potential / eV')
    plt.savefig(plt_fname, dpi=dpi)


def simple_nn(start, elements, end=None, ox_states=None, nn_method=CrystalNN(),
              return_df=False, txt_fname='nn_data.txt'):
    """
    Finds the nearest neighbours for simple structures. Before using on slabs
    make sure the nn_method works with the bulk structure.
    Args:
        start (str): filename of structure, takes all pymatgen-supported formats.
        elements (list): the elements in the structure, order does not matter.
        end (str): filename of structure to analyse, use if comparing initial
        and final structures. The structures must have same constituent atoms
        and number of sites; default=None
        ox_states (list or dict): add oxidation states either by sites
        i.e. [3, 2, 2, 1, -2, -2, -2, -2] or by element i.e. {'Fe': 3, 'O':-2};
        default=None which adds oxidation states by guess
        nn_method (class): the pymatgen.analysis.local_env nearest neighbour
        method; default=CrystalNN()
        return_df (bool): returns the DataFrame; default=False
        txt_fname (str): filename of csv file, default='nn_data.txt'
    Returns
        txt file with atoms and the number of nearest neighbours
    """

    # Instantiate start structure object
    start_struc = Structure.from_file(start)

    # Add atom site labels to the structure
    el_dict = {i : 1 for i in elements}
    site_labels = []
    for site in start_struc:
        symbol = site.specie.symbol
        site_labels.append((symbol,el_dict[symbol]))
        el_dict[symbol] +=1
    start_struc.add_site_property('', site_labels)

    # Adds oxidation states by guess by default or if the oxidation states are
    # antyhing but a list or a dict; max_sites speeds up the by guess method
    if type(ox_states) is dict:
        start_struc.add_oxidation_state_by_element(ox_states)
    elif type(ox_states) is list:
        start_struc.add_oxidation_state_by_site(ox_states)
    else:
        start_struc.add_oxidation_state_by_guess(max_sites=-1)

    # Get bonded start structure
    bonded_start = nn_method.get_bonded_structure(structure=start_struc)

    # Nearest neighbours for just one structure
    if end is None:
        nn_list = []
        for n, site in enumerate(start_struc):
            cn_start = bonded_start.get_coordination_of_site(n)
            label = site_labels[n]
            nn_list.append({'site': n+1,
                           'atom': label,
                           'cn start': cn_start})

        # Make a dataframe from nn_list and return if needed
        df = pd.DataFrame(nn_list)
        df.to_csv(txt_fname, header=True, index=False, sep='\t', mode='w')

        if return_df is True:
            return df

    # Nearest neighbour for two compared structures of the same system
    else:
        end_struc = Structure.from_file(end)

        # Adds oxidation states by guess by default or if the provided oxidation
        # states are antyhing but a list or a dict
        if type(ox_states) is dict:
            end_struc.add_oxidation_state_by_element(ox_states)
        elif type(ox_states) is list:
            end_struc.add_oxidation_state_by_site(ox_states)
        else:
            end_struc.add_oxidation_state_by_guess(max_sites=-1)

        # Get the bonded end structure
        bonded_end = nn_method.get_bonded_structure(end_struc)

        # Get coordination numbers for start and end structures, calculates the
        # end - start difference, collects the site labels and passes it all to
        # a dataframe
        nn_list = []
        for n, site in enumerate(start_struc):
            cn_start = bonded_start.get_coordination_of_site(n)
            cn_end = bonded_end.get_coordination_of_site(n)
            cn_diff = cn_end - cn_start
            label = site_labels[n]
            nn_list.append({'site': n+1,
                           'atom': label,
                           'cn start': cn_start,
                           'cn_end': cn_end,
                           'diff': cn_diff})

        # Make a dataframe from nn_list and return if needed
        df = pd.DataFrame(nn_list)
        df.to_csv(txt_fname, header=True, index=False, sep='\t', mode='w')

        if return_df is True:
            return df


def complex_nn(start, elements, cut_off_dict, end=None, ox_states=None,
               txt_fname='nn_data.txt', return_df=False):
    """
    Finds the nearest neighbours for more complex structures. Uses CutOffDictNN()
    class as the nearest neighbour method. Check validity on bulk structure
    before applying to surface slabs.
    Args:
        start (str): filename of structure to analyse
        elements (list): the elements in the structure, order does not matter
        cut_off_dict (dict): dictionary of bond lengths
        i.e. {('Ag','S'): 3.09, ('La','O'): 2.91, ('La','S'): 3.36,
              ('Ti','O'): 2.35, ('Ti','S'): 2.75, ('Cu','S'): 2.76}
        end (str): filename of structure to analyse, use if comparing initial
        and final structures, the compared structures must have same constituent
        atoms and number of sites; default=None
        ox_states (list or dict): add oxidation states either by sites
        i.e. [3, 2, 2, 1, -2, -2, -2, -2] or by element i.e. {'Fe': 3, 'O':-2}.
        If the structure is decorated with oxidation states, the bond distances
        need to have oxidation states specified. Default=None (no oxidation states)
        txt_fname (str): filename of csv file, default='nn_data.txt'
        return_df (bool): returns the DataFrame; default=False
    Returns
        txt file with atoms and the number of nearest neighbours
    """

     # Instantiate start structure object
    start_struc = Structure.from_file(start)

    # Add atom site labels to the structure
    el_dict = {i : 1 for i in elements}
    site_labels = []
    for site in start_struc:
        symbol = site.specie.symbol
        site_labels.append((symbol,el_dict[symbol]))
        el_dict[symbol] +=1
    start_struc.add_site_property('', site_labels)

    # Adds oxidation states if dict of elements or list of sites are provided,
    # otherwise none are added
    if type(ox_states) is dict:
        start_struc.add_oxidation_state_by_element(ox_states)
    elif type(ox_states) is list:
        start_struc.add_oxidation_state_by_site(ox_states)
    else:
        pass

    # Instantiate the nearest neighbour algorithm
    codnn = CutOffDictNN(cut_off_dict=cut_off_dict)

    # Get the bonded end structure
    bonded_start = codnn.get_bonded_structure(start_struc)

    # Nearest neighbours for just one structure
    if end is None:
        nn_list = []
        for n, site in enumerate(start_struc):
            cn_start = bonded_start.get_coordination_of_site(n)
            label = site_labels[n]
            nn_list.append({'site': n+1,
                           'atom': label,
                           'cn start': cn_start})

        # Make a dataframe from nn_list and return if needed
        df = pd.DataFrame(nn_list)
        df.to_csv(txt_fname, header=True, index=False, sep='\t', mode='w')

        if return_df is True:
            return df

    #nearest neighbour for two compared structures
    else:
        end_struc = Structure.from_file(end)

        # Adds oxidation states if dict of elements or list of sites are provided,
        # otherwise none are added
        if type(ox_states) is dict:
            end_struc.add_oxidation_state_by_element(ox_states)
        elif type(ox_states) is list:
            end_struc.add_oxidation_state_by_site(ox_states)
        else:
            pass

        # Get the bonded end structure
        bonded_end = codnn.get_bonded_structure(end_struc)

        # Get coordination numbers for start and end structures, calculates the
        # end - start difference, collects the site labels and passes it all to
        # a dataframe
        nn_list = []
        for n, site in enumerate(start_struc):
            cn_start = bonded_start.get_coordination_of_site(n)
            cn_end = bonded_end.get_coordination_of_site(n)
            cn_diff = cn_end - cn_start
            label = site_labels[n]
            nn_list.append({'site': n+1,
                           'atom': label,
                           'cn start': cn_start,
                           'cn_end': cn_end,
                           'diff': cn_diff})

        # Make a dataframe from nn_list and return if needed
        df = pd.DataFrame(nn_list)
        df.to_csv(txt_fname, header=True, index=False, sep='\t', mode='w')

        if return_df is True:
            return df

def slab_thickness(start, start_zmax=None, end=None, end_zmax=None):
    """
    Finds the thickness of the slab in c-direction based on the difference
    between the maximum and minimum Cartesian coordinates of the slab.
    Allows slab thickness comparison between the bulk-like and relaxed or
    reconstructed surface slab. Accounts for non-centered slabs and movement of
    atoms to the other side of the slab

    Args:
        start (str): intial structure filename
        start_zmax (float): the Cartesian coordinate of the maximum coordinate
        in c-direction, used if the slab was not centred and some atoms
        moved to the other side of the unit cell; default=None
        end (str): end structure filename, default=None
        end_zmax (float): the cartesian coordinate of the maximum coordinate in
        c-direction, used if the slab was not centred and some atoms moved to
        the other side of the unit cell; default=None

    Returns
        prints slab thickness to terminal
    """
    start_struc = Structure.from_file(start)
    start_struc = start_struc.cart_coords

    # just one slab
    if end is None and start_zmax is None:
        s_xmax, s_ymax, s_zmax = start_struc.max(axis = 0)
        s_xmin, s_ymin, s_zmin = start_struc.min(axis = 0)

        thickness = s_zmax - s_zmin

        print('The slab thickness is {:.3f} Å'.format(thickness))

    # one slab with specified max cartesian coordinate in z
    elif end is None and start_zmax is not None:
        s_xmax, s_ymax, s_zmax = start_struc.max(axis = 0)
        s_xmin, s_ymin, s_zmin = start_struc.min(axis = 0)

        thickness = start_zmax - s_zmin

        print('The slab thickness is {:.3f} Å'.format(thickness))

    # initial and relaxed structures, no spefcified max coords
    elif end_zmax is None:
        end_struc = Structure.from_file(end)
        end_struc = end_struc.cart_coords

        s_xmax, s_ymax, s_zmax = start_struc.max(axis = 0)
        s_xmin, s_ymin, s_zmin = start_struc.min(axis = 0)
        e_xmax, e_ymax, e_zmax = end_struc.max(axis = 0)
        e_xmin, e_ymin, e_zmin = end_struc.min(axis = 0)

        start_thickness = s_zmax - s_zmin
        end_thickness = e_zmax - e_zmin
        difference = start_thickness - end_thickness

        print(('The initial slab thickness is {:.3f} Å. The final slab thickness'
              ' is {:.3f} Å. The difference between the two '
              'is {:.3f} Å').format(start_thickness, end_thickness, difference))

    # intial and relaxed structure, specified max coord
    else:
        end_struc = Structure.from_file(end)
        end_struc = end_struc.cart_coords

        s_xmax, s_ymax, s_zmax = start_struc.max(axis = 0)
        s_xmin, s_ymin, s_zmin = start_struc.min(axis = 0)

        e_xmin, e_ymin, e_zmin = end_struc.min(axis = 0)

        start_thickness = s_zmax - s_zmin
        end_thickness = end_zmax - e_zmin
        difference = start_thickness - end_thickness

        print(('The initial slab thickness is {:.3f} Å. The final slab thickness'
              ' is {:.3f} Å. The difference between the two '
              'is {:.3f} Å').format(start_thickness, end_thickness, difference))
