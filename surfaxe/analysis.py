# Pymatgen
from pymatgen import Structure, Element, Specie
from pymatgen.core.structure import SiteCollection
from pymatgen.core.lattice import Lattice
from pymatgen.analysis.local_env import BrunnerNN_real, BrunnerNN_reciprocal,\
BrunnerNN_relative, CovalentBondNN, Critic2NN, CrystalNN, CutOffDictNN, EconNN,\
JmolNN, MinimumDistanceNN, MinimumOKeeffeNN, MinimumVIRENN, NearNeighbors, VoronoiNN
from pymatgen.io.vasp.outputs import Locpot, VolumetricData

# Misc 
import os
import math
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('once')

# Matplotlib
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
mpl.rcParams['figure.figsize'] = (10.0,8.0)
mpl.rcParams.update({'font.size': 14})

# surfaxe
from surfaxe.generation import oxidation_states
from surfaxe.io import plot_bond_analysis, plot_electrostatic_potential

def cart_displacements(start, end, elements, max_disp=0.1, save_txt=True,
txt_fname='cart_displacements.txt'):
    """
    Produces a text file with all the magnitude of displacements of atoms
    in Cartesian space

    Args:
        start (`str`): Filename of initial structure file in any format 
            supported by pymatgen.
        end (`str`): Filename of final structure file in any format supported
            by pymatgen.
        elements (`list`): List of elements in the structure in any order 
            e.g. ``['Y', 'Ti', 'O', 'S']`` 
        max_disp (`float`, optional): The maximum displacement shown. Defaults 
            to 0.1 Å.
        save_txt (`bool`, optional): Save the displacements to file. Defaults to 
            ``True``.
        txt_fname (`str`, optional): Filename of the csv file. Defaults to 
            ``'cart_displacement.txt'``.

    Returns:
        Displacements of atoms in Cartesian space 

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

    if save_txt: 
        df.to_csv(txt_fname, header=True, index=False, sep='\t', mode='w')
    else: 
        return df

def bond_analysis(structure=None, bonds=None, nn_method=CrystalNN(), 
ox_states=None, save_csv=True, csv_fname='bond_analysis.csv', 
save_plt=True, plt_fname='bond_analysis.png', dpi=300):
    """
    Parses the structure looking for bonds between atoms. Check the validity of
    nearest neighbour method on the bulk structure before using it on slabs.

    Args:
        structure (`str`): filename of structure, takes all pymatgen-supported 
            formats.
        bonds (`list` of `tuples`): List of bonds to compare in any order
            e.g. ``[('Y', 'O'), ('Ti', 'S')]``
        nn_method (`class`, optional): pymatgen.analysis.local_env nearest 
            neighbour method. Defaults to ``CrystalNN()``
        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess. 
              
            Defaults to ``None``. 
        save_csv (`bool`, optional): Makes a csv file with the c coordinate of 
            the first atom and bond length. Defaults to ``True``.
        csv_fname (`str`, optional): Filename of the csv file. Defaults to 
            ``'bond_analysis.csv'``.
        save_plt (`bool`, optional): Make and save the bond analysis plot. 
            Defaults to ``True``. 
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'bond_analysis.png'``.
        dpi (`int`, optional): Dots per inch. Defaults to ``300``. 

    Returns:
        DataFrame with the c coordinate of the first atom and bond length
    """
    # Check all neccessary input parameters are present 
    if not any ([structure, bonds]): 
        raise ValueError('One or more of the required arguments (structure, '
                         'bonds) were not supplied.')

    struc = Structure.from_file(structure)
    struc = oxidation_states(structure=struc, ox_states=ox_states)

    # Iterates through the structure, looking for pairs of bonded atoms. If the
    # sites match, the bond distance is calculated and passed to a dataframe
    bonds_info = []
    for n, pos in enumerate(struc):
        for atom1, atom2 in bonds:
            if pos.specie.symbol is atom1:
                nearest_neighbours = nn_method.get_nn_info(struc, n)
                matched_sites = []
                for d in nearest_neighbours:
                    if d.get('site').specie.symbol is atom2:
                        matched_sites.append(d)
                bond_distances = [struc.get_distance(n,x['site_index']) for x in matched_sites]
                bonds_info.append({
                    '{}_index'.format(atom1): n+1,
                    '{}_c_coord'.format(atom1): pos.c,
                    '{}-{}_bond_distance'.format(atom1,atom2): np.mean(bond_distances)
                })

    df = pd.DataFrame(bonds_info)
    
    # Save plot and csv, or return the DataFrame 
    if save_plt: 
        plot_bond_analysis(bonds, df=df, plt_fname=plt_fname, dpi=dpi)
    if save_csv: 
        df.to_csv(csv_fname, header=True, index=False)
    else: 
        return df


def electrostatic_potential(lattice_vector, locpot='./LOCPOT', axis=2,
save_csv=True, csv_fname='potential.csv', save_plt=True, 
plt_fname='potential.png', dpi=300):
    """
    Reads LOCPOT to get the planar and macroscopic potential in specified 
    direction. 

    Args:
        lattice_vector (`float`): The periodicity of the slab. 
        locpot (`str`, optional): The path to the LOCPOT file. Defaults to 
            ``'./LOCPOT'``
        axis (`int`, optional): The direction in which the potential is 
            investigated; a=0, b=1, c=2. Defaults to `2`. 
        save_csv (`bool`, optional): Makes a csv file with planar and macroscopic 
            potential. Defaults to ``True``.
        csv_fname (`str`, optional): Filename of the csv file. Defaults 
            to ``'potential.csv'``.
        save_plt (`bool`, optional): Make and save the plot of electrostatic 
            potential. Defaults to ``True``. 
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'potential.png'``.
        dpi (`int`, optional): Dots per inch. Defaults to ``300``. 

    Returns:
        DataFrame
    """
    # Read potential and structure data
    lpt = Locpot.from_file(locpot)
    struc = Structure.from_file(locpot)

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

    df = pd.DataFrame(data=planar, columns=['planar'])
    df['macroscopic'] = macroscopic

    # Plot and save the graph, save the csv or return the dataframe
    if save_plt: 
        plot_electrostatic_potential(df=df, plt_fname=plt_fname, dpi=dpi)
    if save_csv: 
        df.to_csv(csv_fname, header=True, index=False)
    else: 
        return df

def simple_nn(start, elements, end=None, ox_states=None, nn_method=CrystalNN(), 
save_csv=True, csv_fname='nn_data.csv'):
    """
    Finds the nearest neighbours for simple structures. Before using on slabs
    make sure the nn_method works with the bulk structure. The required arguments
    are `start` and `elements`. 
    
    Args:
        start (`str`): filename of structure, takes all pymatgen-supported formats.
        elements (`list`): List of elements in the structure in any order 
            e.g. ``['Y', 'Ti', 'O', 'S']`` 
        end (`str`, optional): filename of structure to analyse, use if 
            comparing initial and final structures. The structures must have 
            same constituent atoms and number of sites. Defaults to ``None``. 
        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess. 
              
            Defaults to ``None``. 
        nn_method (`class`, optional): The pymatgen.analysis.local_env nearest 
            neighbour method. Defaults to ``CrystalNN()``.
        save_csv (`bool`, optional): Save to a csv file. Defaults to ``True``.
        csv_fname (`str`, optional): Filename of the csv file. Defaults to 
            ``'nn_data.csv'``
    
    Returns
        DataFrame
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
    
    # Add oxidation states and get bonded structure
    start_struc = oxidation_states(start_struc)
    bonded_start = nn_method.get_bonded_structure(start_struc)

    if end: 
        end_struc = Structure.from_file(end)
        end_struc = oxidation_states(end_struc)
        bonded_end = nn_method.get_bonded_structure(end_struc)
    
    # Iterate through structure, evaluate the coordination number and the 
    # nearest neighbours species for start and end structures, collects the
    # symbol and index of the site (atom) evaluated and its nearest neighbours 
    df_list = []
    for n, site in enumerate(start_struc):
        cn_start = bonded_start.get_coordination_of_site(n)
        coord_start = nn_method.get_nn_info(start_struc, n)
        specie_list = []
        for d in coord_start: 
            spc = d.get('site').specie.symbol 
            specie_list.append(spc)
        specie_list.sort()
        site_nn_start = ' '.join(specie_list)
        label = site_labels[n]

        if end: 
            cn_end = bonded_end.get_coordination_of_site(n)
            coord_end = nn_method.get_nn_info(end_struc, n)
            specie_list = []
            for d in coord_end: 
                spc = d.get('site').specie.symbol 
                specie_list.append(spc)
            specie_list.sort()
            site_nn_end = ' '.join(specie_list)
            df_list.append({'site': n+1, 'atom': label, 'cn start': cn_start,
            'nn_start': site_nn_start, 'cn_end': cn_end, 'nn_end': site_nn_end})

        else: 
            df_list.append({'site_index': n+1, 'site': label,
            'cn_start': cn_start, 'nn_start': site_nn_start})

    # Make a dataframe from df_list 
    df = pd.DataFrame(df_list)

    # Save the csv file or return as dataframe 
    if save_csv: 
        df.to_csv(csv_fname, header=True, index=False)
    else:    
        return df


def complex_nn(start, elements, cut_off_dict, end=None, ox_states=None, 
save_csv=True, csv_fname='nn_data.csv'):
    """
    Finds the nearest neighbours for more complex structures. Uses CutOffDictNN()
    class as the nearest neighbour method. Check validity on bulk structure
    before applying to surface slabs. The required args are `start`, `elements`
    and `cut_off_dict`. 

    Args:
        start (`str`): filename of structure, takes all pymatgen-supported formats.
        elements (`list`): List of elements in the structure in any order 
            e.g. ``['Y', 'Ti', 'O', 'S']`` 
        cut_off_dict (`dict`): Dictionary of bond lengths
            e.g. ``{('Ag','S'): 3.09, ('La', 'O'): 2.91, ('La', 'S'): 3.36,
            ('Ti', 'O'): 2.35, ('Ti', 'S'): 2.75, ('Cu', 'S'): 2.76}``
        end (`str`, optional): filename of structure to analyse, use if 
            comparing initial and final structures. The structures must have 
            same constituent atoms and number of sites. Defaults to ``None``. 
        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: No oxidation states are added - different from other
              functions in ``surfaxe``. 

            Defaults to ``None`` 
        save_csv (`bool`, optional): Save to a csv file. Defaults to ``True``.
        csv_fname (`str`, optional): Filename of the csv file. Defaults to 
            ``'nn_data.csv'``
    
    Returns
        DataFrame
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

    # Add oxidation states 
    if ox_states is not None:
        start_struc = oxidation_states(start_struc, ox_states=ox_states)

    # Instantiate the nearest neighbour algorithm and get bonded structure
    codnn = CutOffDictNN(cut_off_dict=cut_off_dict)
    bonded_start = codnn.get_bonded_structure(start_struc)

    # Instantiate the end structure if provided
    if end: 
        end_struc = Structure.from_file(end)
        end_struc = oxidation_states(end_struc)
        bonded_end = codnn.get_bonded_structure(end_struc)

    # Iterate through structure, evaluate the coordination number and the 
    # nearest neighbours species for start and end structures, collects the
    # symbol and index of the site (atom) evaluated and its nearest neighbours 
    df_list = []
    for n, site in enumerate(start_struc):
        cn_start = bonded_start.get_coordination_of_site(n)
        coord_start = codnn.get_nn_info(start_struc, n)
        specie_list = []
        for d in coord_start: 
            spc = d.get('site').specie.symbol 
            specie_list.append(spc)
        specie_list.sort()
        site_nn_start = ' '.join(specie_list)
        label = site_labels[n]

        if end: 
            cn_end = bonded_end.get_coordination_of_site(n)
            coord_end = codnn.get_nn_info(end_struc, n)
            specie_list = []
            for d in coord_end: 
                spc = d.get('site').specie.symbol 
                specie_list.append(spc)
            specie_list.sort()
            site_nn_end = ' '.join(specie_list)
            df_list.append({'site': n+1, 'atom': label, 'cn start': cn_start,
            'nn_start': site_nn_start, 'cn_end': cn_end, 'nn_end': site_nn_end})

        else: 
            df_list.append({'site_index': n+1, 'site': label,
            'cn_start': cn_start, 'nn_start': site_nn_start})

    # Make a dataframe from df_list 
    df = pd.DataFrame(df_list)
    
    # Save the csv file or return as dataframe 
    if save_csv: 
        df.to_csv(csv_fname, header=True, index=False)
    else:    
        return df

def slab_thickness(start, start_zmax=None, end=None, end_zmax=None):
    """
    Finds the thickness of the slab in c-direction based on the difference
    between the maximum and minimum Cartesian coordinates of the slab.
    Allows slab thickness comparison between the bulk-like and relaxed or
    reconstructed surface slab. Accounts for non-centered slabs and movement of
    atoms to the other side of the slab. The required arg is start. 

    Args:
        start (str): intial structure filename
        start_zmax (float): the Cartesian coordinate of the maximum coordinate
            in c-direction, used if the slab was not centred and some atoms
            moved to the other side of the unit cell; default=None
        end (str): end structure filename, default=None
        end_zmax (float): the cartesian coordinate of the maximum coordinate in
            c-direction, used if the slab was not centred and some atoms moved 
            to the other side of the unit cell; default=None

    Returns
        slab thickness to terminal
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
