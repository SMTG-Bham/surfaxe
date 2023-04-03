# Pymatgen
from pymatgen.core import Structure
from pymatgen.analysis.local_env import CrystalNN, CutOffDictNN
from pymatgen.io.vasp.outputs import Locpot

# Misc
import os
import math
import numpy as np
import pandas as pd
import warnings

# surfaxe
from surfaxe.generation import oxidation_states
from surfaxe.io import plot_bond_analysis, plot_electrostatic_potential, _instantiate_structure

def cart_displacements(start, end, max_disp=0.1, save_txt=True,
txt_fname='cart_displacements.txt'):
    """
    Produces a text file with all the magnitude of displacements of atoms
    in Cartesian space

    Args:
        start (`str`): Filename of initial structure file in any format 
            supported by pymatgen or pymatgen structure object.
        end (`str`): Filename of final structure file in any format supported
            by pymatgen or pymatgen structure object.
        max_disp (`float`, optional): The maximum displacement shown. Defaults 
            to 0.1 Å.
        save_txt (`bool`, optional): Save the displacements to file. Defaults to 
            ``True``.
        txt_fname (`str`, optional): Filename of the csv file. Defaults to 
            ``'cart_displacement.txt'``.

    Returns:
       None (default) or DataFrame of displacements of atoms in Cartesian space 

    """
    # Instantiate the structures 
    start_struc = _instantiate_structure(start)
    end_struc = _instantiate_structure(end)

    # Add the site labels to the structure
    els = ''.join([i for i in start_struc.formula if not i.isdigit()]).split(' ')
    el_dict = {i : 1 for i in els}
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
            disp_list.append({
                'site': n+1,
                'atom': label,
                # this makes the displacements round to the same number of 
                # decimal places as max displacement, for presentation 
                'displacement': round(d, int(format(max_disp, 'E')[-1])) 
                             })
    # Save as txt file
    df = pd.DataFrame(disp_list)

    if save_txt: 
        df.to_csv(txt_fname, header=True, index=False, sep='\t', mode='w')
    else: 
        return df

def bond_analysis(structure, bond, nn_method=CrystalNN(), ox_states=None, 
save_csv=True, csv_fname='bond_analysis.csv', save_plt=False, 
plt_fname='bond_analysis.png', **kwargs):
    """
    Parses the structure looking for bonds between atoms. Check the validity of
    the nearest neighbour method on the bulk structure before using it on slabs.

    Args:
        structure (`str`): filename of structure, takes all pymatgen-supported 
            formats, including pmg structure object
        bond (`list`): Bond to analyse e.g. ``['Y', 'O']``
        nn_method (`class`, optional): The coordination number prediction 
            algorithm used. Because the ``nn_method`` is a class, the class 
            needs to be imported from ``pymatgen.analysis.local_env`` before it 
            can be instantiated here. Defaults to ``CrystalNN()``.
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
            Defaults to ``False``. 
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'bond_analysis.png'``. 
        kwargs (`dict`, optional): Additional keyword arguments to pass to 
            plot_bond_analysis(). Defaults to ``{}``.

    Returns:
        DataFrame with the c coordinate of the first atom and bond length
    """
    struc = _instantiate_structure(structure)
    struc = oxidation_states(structure=struc, ox_states=ox_states)

    if len(bond) > 2: 
        warnings.warn('Bond with more than two elements supplied. '
        'Only the first two elements will be treated as a bond.')

    # Iterates through the structure, looking for pairs of bonded atoms. If the
    # sites match, the bond distance is calculated and passed to a dataframe
    bonds_info = []
    for n, pos in enumerate(struc):
        if pos.specie.symbol == bond[0]:
            nearest_neighbours = nn_method.get_nn_info(struc, n)
            matched_sites = []
            for d in nearest_neighbours:
                if d.get('site').specie.symbol == bond[1]:
                    matched_sites.append(d)
            bond_distances = [
                struc.get_distance(n,x['site_index']) for x in matched_sites
            ]
            bonds_info.append({
                '{}_index'.format(bond[0]): n+1,
                '{}_c_coord'.format(bond[0]): pos.c,
                '{}-{}_bond_distance'.format(bond[0],bond[1]): np.mean(bond_distances)
            })

    df = pd.DataFrame(bonds_info)
    
    # Save plot and csv, or return the DataFrame 
    if save_plt: 
        plot_bond_analysis(bond, df=df, plt_fname=plt_fname, **kwargs)
    if save_csv: 
        if not csv_fname.endswith('.csv'):
            csv_fname += '.csv'
        df.to_csv(csv_fname, header=True, index=False)
    else: 
        return df


def electrostatic_potential(locpot='./LOCPOT', prim_to_conv=1,
axis='c', save_csv=True, csv_fname='potential.csv', save_plt=True, 
plt_fname='potential.png', lattice_vector=None, **kwargs):
    """
    Reads LOCPOT to get the planar and macroscopic potential in a, b or c direction 

    Args:
        locpot (`str`, optional): The path to the LOCPOT file. Defaults to 
            ``'./LOCPOT'``
        prim_to_conv (`int`, optional): The number of primitive cells in the 
            conventional cell. Defaults to ``1``.
        axis (`str`, optional): Axis along which the potential is calculated. 
            Takes a,b,c or x,y,z. 
        save_csv (`bool`, optional): Saves to csv. Defaults to ``True``.
        csv_fname (`str`, optional): Filename of the csv file. Defaults 
            to ``'potential.csv'``.
        save_plt (`bool`, optional): Make and save the plot of electrostatic 
            potential. Defaults to ``True``. 
        plt_fname (`str`, optional): Filename of the plot. Defaults to 
            ``'potential.png'``.
        lattice_vector (`float`, optional): Manually set the periodicity of the slab

    Returns:
        DataFrame
    """
    # set up axis 
    if axis in ['a', 'x']: 
        ax = 0 
    elif axis in ['b', 'y']: 
        ax = 1
    elif axis in ['c', 'z']: 
        ax = 2
    else: 
        raise ValueError('axis can only be set to a,b,c or x,y,z')

    # Read potential and structure data
    if os.path.exists(locpot):
        lpt = Locpot.from_file(locpot)
        struc = lpt.structure
    elif os.path.exists(locpot + ".gz"):
        lpt = Locpot.from_file(locpot + ".gz")
        struc = lpt.structure
    else:
        raise FileNotFoundError(
            f"""No LOCPOT(.gz) found at {locpot}(.gz)""")

    # Planar potential
    planar = lpt.get_average_along_axis(ax)
    df = pd.DataFrame(data=planar, columns=['planar']) 
    
    # Calculate macroscopic potential
    if lattice_vector is None:
        # Calculate lattice vector
        arr = np.array([i.coords for i in struc.sites])
        comp, factor = struc.composition.get_reduced_composition_and_factor()
        argmin = arr[:, ax].argmin() 
        specie_min = str(struc[argmin].specie)
        argmax = int(comp.as_dict()[specie_min] * prim_to_conv + argmin)
        lattice_vector = arr[:, ax][argmax] - arr[:, ax][argmin]

    # Divide lattice parameter by no. of grid points in the direction
    resolution = struc.lattice.abc[ax]/lpt.dim[ax]

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
    df['macroscopic'] = macroscopic

    # Get gradient of the plot - this is used for convergence testing, to make 
    # sure the potential is actually flat
    df['gradient'] = np.gradient(df['planar'])

    # Plot and save the graph, save the csv or return the dataframe
    if save_plt: 
        plot_electrostatic_potential(df=df, plt_fname=plt_fname, **kwargs)
    if save_csv: 
        if not csv_fname.endswith('.csv'):
            csv_fname += '.csv'
        df.to_csv(csv_fname, header=True, index=False)
    else: 
        return df

def surface_dipole(filename, **kwargs): 
    """
    Calculates surface dipole for a slab. Useful for band alignments. 

    Args:
        filename (`str`): The path to the LOCPOT or the csv file with 
            macroscopic potential
        kwargs: Keyword arguments for ``electrostatic_potential``
    Returns:
        float: The surface dipole in eV
    """
    if 'LOCPOT' in filename: 
        pt = electrostatic_potential(filename, **kwargs) 
    elif filename.endswith('csv'): 
        pt = pd.read_csv(filename) 
        if 'macroscopic' not in pt.columns: 
            raise ValueError('csv should contain macroscopic potential') 
    else: 
        raise ValueError('filename should be a LOCPOT or a csv file')
    
    # Get the surface dipole 
    dipole = pt['macroscopic'].max() - pt['macroscopic'].iloc[int(len(pt)/2)]
    return round(dipole, 3)

def simple_nn(start, end=None, ox_states=None, nn_method=CrystalNN(), 
save_csv=True, csv_fname='nn_data.csv'):
    """
    Finds the nearest neighbours for simple structures. Before using on slabs
    make sure the nn_method works with the bulk structure. 
    
    The ``site_index`` in the produced DataFrame or csv file is one-indexed and 
    represents the atom index in the structure. 

    Args:
        start (`str`): Filename of structure file in any format supported by 
            pymatgen
        end (`str`, optional): Filename of structure file in any format 
            supported by pymatgen. Use if comparing initial and final structures. 
            The structures must have same constituent atoms and number of sites. 
            Defaults to ``None``. 
        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess. 
              
            Defaults to ``None``. 
        nn_method (`class`, optional): The coordination number prediction 
            algorithm used. Because the ``nn_method`` is a class, the class 
            needs to be imported from pymatgen.analysis.local_env before it 
            can be instantiated here. Defaults to ``CrystalNN()``.
        save_csv (`bool`, optional): Save to a csv file. Defaults to ``True``.
        csv_fname (`str`, optional): Filename of the csv file. Defaults to 
            ``'nn_data.csv'``
    
    Returns
        None (default) or DataFrame containing coordination data 
    """
    # Instantiate start structure object
    start_struc = _instantiate_structure(start)

    # Add atom site labels to the structure
    els = ''.join([i for i in start_struc.formula if not i.isdigit()]).split(' ')
    el_dict = {i : 1 for i in els}
    site_labels = []
    for site in start_struc:
        symbol = site.specie.symbol
        site_labels.append((symbol,el_dict[symbol]))
        el_dict[symbol] +=1
    start_struc.add_site_property('', site_labels)
    
    # Add oxidation states and get bonded structure
    start_struc = oxidation_states(start_struc, ox_states)
    bonded_start = nn_method.get_bonded_structure(start_struc)

    if end: 
        end_struc = _instantiate_structure(end)
        end_struc = oxidation_states(end_struc, ox_states)
        bonded_end = nn_method.get_bonded_structure(end_struc)
    
    # Iterate through structure, evaluate the coordination number and the 
    # nearest neighbours specie for start and end structures, collects the
    # symbol and index of the site (atom) evaluated and its nearest neighbours 
    df_list = []
    for n, site in enumerate(start_struc):
        cn_start = bonded_start.get_coordination_of_site(n)
        coord_start = bonded_start.get_connected_sites(n)
        specie_list = []
        for d in coord_start: 
            spc = d.site.specie.symbol 
            specie_list.append(spc)
        specie_list.sort()
        site_nn_start = ' '.join(specie_list)
        label = site_labels[n]

        if end: 
            cn_end = bonded_end.get_coordination_of_site(n)
            coord_end = bonded_end.get_connected_sites(n)
            specie_list = []
            for d in coord_end: 
                spc = d.site.specie.symbol 
                specie_list.append(spc)
            specie_list.sort()
            site_nn_end = ' '.join(specie_list)
            df_list.append({'site': n+1, 'atom': label, 'cn_start': cn_start,
            'nn_start': site_nn_start, 'cn_end': cn_end, 'nn_end': site_nn_end})

        else: 
            df_list.append({'site_index': n+1, 'site': label,
            'cn_start': cn_start, 'nn_start': site_nn_start})

    # Make a dataframe from df_list 
    df = pd.DataFrame(df_list)

    # Save the csv file or return as dataframe 
    if save_csv: 
        if not csv_fname.endswith('.csv'):
            csv_fname += '.csv'
        df.to_csv(csv_fname, header=True, index=False)
    else:    
        return df


def complex_nn(start,  cut_off_dict, end=None, ox_states=None, 
save_csv=True, csv_fname='nn_data.csv'):
    """
    Finds the nearest neighbours for more complex structures. Uses CutOffDictNN()
    class as the nearest neighbour method. Check validity on bulk structure
    before applying to surface slabs. 

    The ``site_index`` in the produced DataFrame or csv file is one-indexed and 
    represents the atom index in the structure. 

    Args:
        start (`str`): filename of structure, takes all pymatgen-supported formats.
        cut_off_dict (`dict`): Dictionary of bond lengths. The bonds should be 
            specified with the oxidation states\n
            e.g. ``{('Bi3+', 'O2-'): 2.46, ('V5+', 'O2-'): 1.73}``
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
            
            * if ``None``:  The oxidation states are added by guess.  

            Defaults to ``None`` 
        save_csv (`bool`, optional): Save to a csv file. Defaults to ``True``.
        csv_fname (`str`, optional): Filename of the csv file. Defaults to 
            ``'nn_data.csv'``
    
    Returns
        None (default) or DataFrame containing coordination data.
    """
    # Instantiate start structure object
    start_struc = Structure.from_file(start)

    # Add atom site labels to the structure
    els = ''.join([i for i in start_struc.formula if not i.isdigit()]).split(' ')
    el_dict = {i : 1 for i in els}
    site_labels = []
    for site in start_struc:
        symbol = site.specie.symbol
        site_labels.append((symbol,el_dict[symbol]))
        el_dict[symbol] +=1
    start_struc.add_site_property('', site_labels)

    # Add oxidation states 
    start_struc = oxidation_states(start_struc, ox_states=ox_states)

    # Instantiate the nearest neighbour algorithm and get bonded structure
    codnn = CutOffDictNN(cut_off_dict=cut_off_dict)
    bonded_start = codnn.get_bonded_structure(start_struc)

    # Instantiate the end structure if provided
    if end: 
        end_struc = Structure.from_file(end)
        end_struc = oxidation_states(end_struc, ox_states=ox_states)
        bonded_end = codnn.get_bonded_structure(end_struc)

    # Iterate through structure, evaluate the coordination number and the 
    # nearest neighbours specie for start and end structures, collects the
    # symbol and index of the site (atom) evaluated and its nearest neighbours 
    df_list = []
    for n, site in enumerate(start_struc):
        cn_start = bonded_start.get_coordination_of_site(n)
        coord_start = bonded_start.get_connected_sites(n)
        specie_list = []
        for d in coord_start: 
            spc = d.site.specie.symbol 
            specie_list.append(spc)
        specie_list.sort()
        site_nn_start = ' '.join(specie_list)
        label = site_labels[n]

        if end: 
            cn_end = bonded_end.get_coordination_of_site(n)
            coord_end = bonded_end.get_connected_sites(n)
            specie_list = []
            for d in coord_end: 
                spc = d.site.specie.symbol 
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
        if not csv_fname.endswith('.csv'):
            csv_fname += '.csv'
        df.to_csv(csv_fname, header=True, index=False)
    else:    
        return df
