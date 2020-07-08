
from pymatgen.core.trajectory import Trajectory
from pymatgen import Structure, Element, Specie
from pymatgen.core.structure import SiteCollection
from pymatgen.analysis.local_env import CrystalNN, BrunnerNN_reciprocal, CutOffDictNN, VoronoiNN, MinimumDistanceNN
import os
import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import warnings
warnings.filterwarnings('once')
import matplotlib as mpl
### Matplotlib style ###
mpl.rcParams['figure.figsize'] = (7.0,6.0)


def cart_displacements(start, end, elements, max_disp=0.1):
    """
    gives a text? file with all the cartesian displacements of atoms

    Args:
        start (str): initial structure
        end (str): final structure
        elements (list): list of elements in the structure e.g. ['Y', 'Ti', 'O', 'S']
        max_disp (float): maximum displacement shown, default 0.1 Å
    Returns:
        text file w displacements

    """

    start_struc = Structure.from_file(start)
    end_struc = Structure.from_file(end)

    el_dict = {i : 1 for i in elements}
    site_labels = []

    for site in start_struc:
        symbol = site.specie.symbol
        site_labels.append((symbol,el_dict[symbol]))
        el_dict[symbol] +=1
    start_struc.add_site_property('', site_labels)

    start_struc = start_struc.cart_coords
    end_struc = end_struc.cart_coords

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

    df = pd.DataFrame(disp_list)
    df.to_csv('cart_displacements.txt', header = True, index = False, sep='\t', mode='w')

#raise warning if start, end don't have the same number of sites; raise warning if all disp smaller than the max_disp;

def slab_thickness(start, end=None, end_zmax=None):
    """
    finds the thickness of the slab in c-direction based on the difference between the maximum and minimum cartesian coordinates of the slab
    allows for slab thickness comparison between the bulk-like and relaxed/reconstructed surface slab by setting end to final structure filename
    Args:
        start (str?): intial structure filename
        end (str): end structure filename, default: None
        end_zmax (float): the cartesian coordinate of the maximum coordinate in c-direction, used if the slab was not centred and some atoms moved to the other side of the unit cell; default: None
    Returns
        slab thickness
    """
    start_struc = Structure.from_file(start)
    start_struc = start_struc.cart_coords

    if end is None:
        s_xmax, s_ymax, s_zmax = start_struc.max(axis = 0)
        s_xmin, s_ymin, s_zmin = start_struc.min(axis = 0)

        thickness = s_zmax - s_zmin

        print('the slab thickness is {:.3f}'.format(thickness))

    elif end_zmax is None:
        end_struc = Structure.from_file(end)
        end_struc = contcar.cart_coords

        s_xmax, s_ymax, s_zmax = start_struc.max(axis = 0)
        s_xmin, s_ymin, s_zmin = start_struc.min(axis = 0)
        e_xmax, e_ymax, e_zmax = end_struc.max(axis = 0)
        e_xmin, e_ymin, e_zmin = end_struc.min(axis = 0)

        start_thickness = s_zmax - s_zmin
        end_thickness = e_zmax - e_zmin
        difference = start_thickness - end_thickness

        print('the initial slab thickness is {:.3f}. the final slab thickness is {:.3f} and the difference between the two is {:.3f}'.format(start_thickness, end_thickness, difference))

    else:
        end_struc = Structure.from_file(end)
        end_struc = contcar.cart_coords

        s_xmax, s_ymax, s_zmax = start_struc.max(axis = 0)
        s_xmin, s_ymin, s_zmin = start_struc.min(axis = 0)

        e_xmin, e_ymin, e_zmin = end_struc.min(axis = 0)

        start_thickness = s_zmax - s_zmin
        end_thickness = end_zmax - e_zmin
        difference = start_thickness - end_thickness

        print('the initial slab thickness is {:.3f}. the final slab thickness is {:.3f} and the difference between the two is {:.3f}'.format(start_thickness, end_thickness, difference))



def bond_analysis(struc, atoms, nn_method=VoronoiNN(tol=0.1, cutoff=10)):
    """
    parses the structure looking for bonds between atoms
    Args:
        struc: filename or structure file
        atoms (list of tuples): list of bonds to compare e.g. [('Y', 'O'), ('Ti', 'S')]
        nn_method (class?): nearest neighbour method to use to find the coordination environment
    Returns:
        csv file
    """

    bonds_info = []
    struc.add_oxidation_state_by_guess()

    for n, pos in enumerate(struc):
        for atom1, atom2 in atoms:
            if pos.specie.symbol == atom1:
                nearest_neighbours =  nn_method.get_nn_info(struc, n)
                matched_sites = []
                for d in nearest_neighbours:
                    if d.get('site').specie.symbol == atom2:
                        matched_sites.append(d)
                bond_distances = [struc.get_distance(n,x['site_index']) for x in matched_sites]
                bonds_info.append({'{}_index'.format(atom1): n+1,
                                   '{}_c_coord'.format(atom1): pos.c,
                                   '{}-{}_bond_distance'.format(atom1,atom2): np.mean(bond_distances)})
    df = pd.DataFrame(bonds_info)
    df.to_csv('bond_analysis_data.csv', index=False)

    return df


def plot_bond_analysis(atoms):
    """
    plots the bond distance graph thingie luisa introduced to the mix from the csv file generated with bond_analysis
    Args:
        atoms (list of tuples) in the same order as atoms in bond_analysis
    Returns:
        plot of bond distance change in c direction
    """

    df = pd.read_csv('bond_analysis_data.csv')
    colors = plt.rcParams["axes.prop_cycle"]()

    legend_list = []
    fig, axs = plt.subplots(len(atoms))

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
    plt.show()


#this is still very much a work in progress, i cannot figure out how to make the fkcin nearest neighbour methods work and i am too scared to ask alex again lol
def nearest_neighbour(start, elements, end=None, nn_method=MinimumDistanceNN()):
    """

    Args:
        start (str):
        elements (list):
        end (str): default None
        nn_method :
    """

    start_struc = Structure.from_file(start)

    el_dict = {i : 1 for i in elements}
    site_labels = []

    for site in start_struc:
        symbol = site.specie.symbol
        site_labels.append((symbol,el_dict[symbol]))
        el_dict[symbol] +=1
    start_struc.add_site_property('', site_labels)

    #start_struc.merge_sites(tol=0.001)

    start_struc.add_oxidation_state_by_guess()


    if end is None:
        nn_list = []
        for n, site in enumerate(start_struc):
            cn_start = nn_method.get_cn(start_struc, n)
            label = site_labels[n]
            nn_list.append({'site': n+1,
                           'atom': label,
                           'cn start': cn_start})
        df = pd.DataFrame(nn_list)
        df.to_csv('nn_data.txt', header = True, index = False, sep='\t', mode='w')

    else:
        end_struc = Structure.from_file(end)
        end_struc.add_oxidation_state_by_guess()

        nn_list = []
        for n, waa in enumerate(start_struc):
            cn_start = nn_method.get_cn(start_struc, n)
            cn_end = nn_method.get_cn(end_struc, n)
            cn_diff = cn_end - cn_start
            label = site_labels[n]
            nn_list.append({'site': n+1,
                           'atom': label,
                           'cn start': cn_start,
                           'cn_end': cn_end,
                           'diff': cn_diff})
        df = pd.DataFrame(nn_list)
        df.to_csv('nn_data.txt', header = True, index = False, sep='\t', mode='w')


# everything below this is my og very dodgy nearest neighbour analysis tool thing it works but i've been trying to make it work for a general case
poscar = Structure.from_file('POSCAR_LTA_010')
poscar.add_oxidation_state_by_guess

contcar = Structure.from_file('CONTCAR_LTA_010')
contcar.add_oxidation_state_by_guess

cod = {
    ('Ag','S'): 3.09,
    ('La','O'): 2.91,
    ('La','S'): 3.36,
    ('Ti','O'): 2.35,
    ('Ti','S'): 2.75,
    ('Cu','S'): 2.76

}

codnn = CutOffDictNN(cut_off_dict=cod)

el_dict = {'La':1 ,'Ti':1 ,'Ag':1 ,'S':1 ,'O':1}
site_labels = []

for site in poscar:
    symbol = site.specie.symbol
    site_labels.append((symbol,el_dict[symbol]))
    el_dict[symbol] +=1
poscar.add_site_property('', site_labels)

for n, site in enumerate(poscar):
    cn_pos = codnn.get_cn(poscar,n)
    cn_cont = codnn.get_cn(contcar,n)
    print(n+1, site.properties.get(''), cn_pos, cn_cont, cn_cont-cn_pos)
