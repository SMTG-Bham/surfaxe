# pymatgen
from pymatgen.core.surface import SlabGenerator, get_symmetrically_distinct_miller_indices
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# misc
import warnings
import itertools
import functools
import multiprocessing
import math
import numpy as np 
import json
from copy import deepcopy

# surfaxe
from surfaxe.io import slabs_to_file, _custom_formatwarning, _instantiate_structure

def generate_slabs(structure, hkl, thicknesses, vacuums, save_slabs=True, 
save_metadata=True, json_fname=None, make_fols=False, make_input_files=False, 
max_size=500, center_slab=True, ox_states=None, is_symmetric=True, 
layers_to_relax = None, fmt='poscar', name='POSCAR', config_dict=None, 
user_incar_settings=None, user_kpoints_settings=None, user_potcar_settings=None, 
parallelise=True, processes=None, **kwargs): 
    """
    Generates all unique slabs for a specified Miller indices or up to a maximum 
    Miller index with minimum slab and vacuum thicknesses. It includes all 
    combinations for multiple zero dipole symmetric terminations for 
    the same Miller index. 
    
    The function returns None by default and generates either: 

    (i) POSCAR_hkl_slab_vac_index.vasp (default) 
    (ii) hkl/slab_vac_index folders with structure files
    (iii) hkl/slab_vac_index with all VASP input files 
    
    Or if `save_slabs=False` a list of dicts of all unique slabs is returned. 
    
    Args:
        structure (`str` or pmg Structure obj): Filename of structure file in 
            any format supported by pymatgen or pymatgen structure object. 
        hkl (`tuple`, `list` or `int`): Miller index as tuple, a list of Miller 
            indices or a maximum index up to which the search should be 
            performed. E.g. if searching for slabs up to (2,2,2) ``hkl=2``
        thicknesses (`list`): The minimum size of the slab in Angstroms. 
        vacuums (`list`): The minimum size of the vacuum in Angstroms.
        save_slabs (`bool`, optional): Whether to save the slabs to file. 
            Defaults to ``True``.
        save_metadata (`bool`, optional): Whether to save the slabs' metadata to 
            file. Saves the entire slab object to a json file 
            Defaults to ``True``. 
        json_fname (`str`, optional): Filename of json metadata file. Defaults to
            bulk_formula_metadata.json
        make_fols (`bool`, optional): Makes folders for each termination 
            and slab/vacuum thickness combinations containing structure files. 
            
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

            Defaults to ``False``.    
        make_input_files (`bool`, optional): Makes INCAR, POTCAR and 
            KPOINTS files in each folder. If ``make_input_files`` is ``True`` 
            but ``make_files`` or ``save_slabs`` is ``False``, files will be 
            saved to folders regardless. This only works with VASP input files, 
            other formats are not yet supported. Defaults to ``False``. 
        max_size (`int`, optional): The maximum number of atoms in the slab 
            specified to raise warning about slab size. Even if the warning is 
            raised, it still outputs the slabs regardless. Defaults to ``500``. 
        center_slab (`bool`, optional): The position of the slab in the 
            simulation cell. 
            
            * ``True``: the slab is centered with equal amounts of 
              vacuum above and below.

            * ``False``: the slab is at the bottom of the simulation cell with
              all of the vacuum on top of it. 

            Defaults to True. 

        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the bulk structure. Different types of oxidation states specified 
            will result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess. 
              
            Defaults to ``None``. 
        is_symmetric (`bool`, optional): Whether the slabs cleaved should 
            have inversion symmetry. If bulk is non-centrosymmetric, 
            ``is_symmetric`` needs to be ``False`` - the function will return no
            slabs as it looks for inversion symmetry. Take care checking the 
            slabs for mirror plane symmetry before just using them. Defaults to 
            ``True``. 
        layers_to_relax (`int`, optional): Specifies the number of layers at the 
            top and bottom of the slab that should be relaxed, keeps the centre
            constrained using selective dynamics. NB only works for VASP files 
        fmt (`str`, optional): The format of the output structure files. Options 
            include 'cif', 'poscar', 'cssr', 'json', not case sensitive. 
            Defaults to 'poscar'. 
        name (`str`, optional): The name of the surface slab structure file 
            created. Case sensitive. Defaults to 'POSCAR'
        config_dict (`dict` or `str`, optional): Specifies the dictionary used 
            for the generation of the input files. Suppports already loaded 
            dictionaires, yaml and json files. Surfaxe-supplied dictionaries 
            are PBE (``pe``), PBEsol (``ps``) and HSE06 (``hse06``) for single 
            shot calculations and PBE (``pe_relax``) and PBEsol (``ps_relax``) 
            for relaxations. Not case sensitive. Defaults to PBEsol (``ps``). 
        user_incar_settings (`dict`, optional): Overrides the default INCAR 
            parameter settings. Defaults to ``None``.
        user_kpoints_settings (`dict` or Kpoints object, optional): 
            Overrides the default kpoints settings. If it is supplied  
            as `dict`, it should be as ``{'reciprocal_density': 100}``. Defaults 
            to ``None``.
        user_potcar_settings (`dict`, optional): Overrides the default POTCAR 
            settings. Defaults to ``None``.
        parallelise (`bool`, optional): Use multiprocessing to generate
            slabs. Defaults to ``True``. 
        processes (`int`, optional): Number of cpu processes to use for     
            multiprocessing. Defauts to ``None``, which is max-1 available.

    Returns:
        None (default) 
        or unique_slabs (list of dicts) 
    """

    # Set up additional arguments for multiprocessing and saving slabs
    mp_kwargs = {'in_unit_planes': False, 'primitive': True, 
    'max_normal_search': None, 'reorient_lattice': True, 'lll_reduce': True, 
    'ftol': 0.1, 'tol': 0.1, 'max_broken_bonds': 0, 'symmetrize': False, 
    'repair': False, 'bonds': None}
    mp_kwargs.update(
        (k, kwargs[k]) for k in mp_kwargs.keys() & kwargs.keys()
    )

    save_slabs_kwargs = {'user_incar_settings': None, 
    'user_kpoints_settings': None, 'user_potcar_settings': None, 
    'constrain_total_magmom': False, 'sort_structure': True, 
    'potcar_functional': None, 'user_potcar_functional': None, 
    'force_gamma': False, 'reduce_structure': None, 'vdw': None, 
    'use_structure_charge': False, 'standardize': False, 'sym_prec': 0.1, 
    'international_monoclinic': True}
    save_slabs_kwargs.update(
        (k, kwargs[k]) for k in save_slabs_kwargs.keys() & kwargs.keys() 
    )
    save_slabs_kwargs.update({'user_incar_settings': user_incar_settings, 
        'user_kpoints_settings': user_kpoints_settings, 
        'user_potcar_settings': user_potcar_settings})
    
    # Set up multiprocessing
    if processes == None:
        processes = multiprocessing.cpu_count() - 1

    # Import bulk relaxed structure, add oxidation states for slab dipole
    # calculations
    struc = _instantiate_structure(structure)

    # Check structure is conventional standard and warn if not 
    sga = SpacegroupAnalyzer(struc)
    cs = sga.get_conventional_standard_structure()
    if cs.lattice != struc.lattice: 
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Lattice of the structure provided does not match the '
        'conventional standard structure. Miller indices are lattice dependent,'
        ' make sure you are using the correct bulk structure ')
    
    # Check if oxidation states were added to the bulk already 
    struc = oxidation_states(struc, ox_states)
    
    # Check if hkl provided as tuple or int, find all available hkl if
    # provided as int; make into a list to iterate over
    if type(hkl) == tuple: 
        miller = [hkl]
    elif type(hkl) == int: 
        miller = get_symmetrically_distinct_miller_indices(struc, hkl)
    elif type(hkl) == list and all(isinstance(x, tuple) for x in hkl): 
        miller = hkl 
    else: 
        raise TypeError('Miller index should be supplied as tuple, int or list '
        'of tuples')
    
    # create all combinations of hkl, slab and vacuum thicknesses
    combos = itertools.product(miller, thicknesses, vacuums)

    # Check if bulk structure is noncentrosymmetric if is_symmetric=True, 
    # change to False if not to make sure slabs are produced, issues warning 
    if is_symmetric: 
        sg = SpacegroupAnalyzer(struc)
        if not sg.is_laue(): 
            is_symmetric = False
            warnings.formatwarning = _custom_formatwarning
            warnings.warn(('Inversion symmetry was not found in the bulk '
            'structure, slabs produced will be non-centrosymmetric'))
    
    # Check if multiple cores are available, then iterate through the slab and 
    # vacuum thicknesses and get all non polar symmetric slabs  
    if multiprocessing.cpu_count() > 1 and parallelise==True:
        with multiprocessing.Pool(processes) as pool:
            nested_provisional = pool.starmap(
                    functools.partial(_mp_generate_slabs, struc,
                    is_symmetric=is_symmetric, center_slab=center_slab,
                    **mp_kwargs), combos)

        provisional = list(itertools.chain.from_iterable(nested_provisional)) 

    else: 
        # Set up kwargs again 
        SG_kwargs = {k: mp_kwargs[k] for k in ['in_unit_planes', 'primitive', 'max_normal_search', 'reorient_lattice', 'lll_reduce']}
        gs_kwargs = {k: mp_kwargs[k] for k in ['ftol', 'tol', 'max_broken_bonds', 
        'symmetrize', 'repair', 'bonds']}

        provisional = []
        for hkl, thickness, vacuum in combos:
            slabgen = SlabGenerator(struc, hkl, thickness, vacuum,
                                    center_slab=center_slab,  
                                    **SG_kwargs) 

            # Get the number of layers in the slab                        
            h = slabgen._proj_height
            p = round(h/slabgen.parent.lattice.d_hkl(slabgen.miller_index), 8)
            if slabgen.in_unit_planes:
                nlayers_slab = int(math.ceil(slabgen.min_slab_size / p))
            else: 
                nlayers_slab = int(math.ceil(slabgen.min_slab_size / h))

            slabs = slabgen.get_slabs(**gs_kwargs)
            for i, slab in enumerate(slabs):
                # Get all the zero-dipole slabs with inversion symmetry
                if is_symmetric: 
                    if slab.is_symmetric() and not slab.is_polar():
                        provisional.append({
                            'hkl': ''.join(map(str, slab.miller_index)),
                            'slab_thickness': thickness,
                            'slab_layers': nlayers_slab,
                            'vac_thickness': vacuum,
                            'slab_index': i,
                            'slab': slab})
                
                # Get all the zero-dipole slabs wihtout inversion symmetry
                else: 
                    if not slab.is_polar():
                        provisional.append({
                            'hkl': ''.join(map(str, slab.miller_index)),
                            'slab_thickness': thickness,
                            'slab_layers': nlayers_slab,
                            'vac_thickness': vacuum,
                            'slab_index': i,
                            'slab': slab})
                  
    # Iterate though provisional slabs to extract the unique slabs
    unique_list_of_dicts, repeat, large = _filter_slabs(provisional, max_size)

    if layers_to_relax is not None and fmt.lower() == 'poscar': 
        unique_list_of_dicts, small = _get_selective_dynamics(struc, 
        unique_list_of_dicts, layers_to_relax)        

        if small: 
            warnings.formatwarning = _custom_formatwarning
            warnings.warn('Some slabs were too thin to fix the centre of the slab.'
            ' Slabs with no selective dynamics applied are: ' + 
            ', '.join(map(str, small)))

    # Warnings for too large, too small, repeated and no slabs; 
    if repeat:
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Not all combinations of hkl or slab/vac thicknesses '
        'were generated because of repeat structures. '
        'The repeat slabs are: ' + ', '.join(map(str, repeat)))

    if large:
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Some generated slabs exceed the max size specified.'
        ' Slabs that exceed the max size are: ' + ', '.join(map(str, large)))
    
    if len(unique_list_of_dicts) == 0: 
        raise ValueError('No zero dipole (Tasker I or II) slabs found for specified Miller index')

    # Save the metadata or slabs to file or return the list of dicts 
    if save_metadata: 
        bulk_name = struc.composition.reduced_formula
        if json_fname is None: 
            json_fname = '{}_metadata.json'.format(bulk_name)
        unique_list_of_dicts_copy = deepcopy(unique_list_of_dicts)
        for i in unique_list_of_dicts_copy: 
            i['slab'] = i['slab'].as_dict() 
        with open(json_fname, 'w') as f: 
            json.dump(unique_list_of_dicts_copy, f)

    if save_slabs: 
        slabs_to_file(list_of_slabs=unique_list_of_dicts, structure=struc, 
        make_fols=make_fols, make_input_files=make_input_files, 
        config_dict=config_dict, fmt=fmt, name=name, **save_slabs_kwargs)
    
    else: 
        return unique_list_of_dicts

def oxidation_states(structure, ox_states=None):
    ''' 
    Adds oxidation states to the structure object if not already present
    
    Args: 
        structure (`obj`): Pymatgen structure object
        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess. 

    Returns: 
        Structure decorated with oxidation states 
    ''' 
    try: 
        structure.species[0].oxi_state
    
    except AttributeError: 
        if type(ox_states) is dict:
            structure.add_oxidation_state_by_element(ox_states)
        elif type(ox_states) is list:
            structure.add_oxidation_state_by_site(ox_states)
        else:
            structure.add_oxidation_state_by_guess(max_sites=-1)

    return structure

def _filter_slabs(provisional, max_size): 
    """
    Filters the repeat slabs from the list of all the zero dipole slabs. 
    Creates lists of large and repeat slabs if any are present for the warnings

    Args: 
        provisional (`list`): All zero dipole slabs generated with SlabGenerator
        max_size (`int`): The maximum number of atoms in the slab 
            specified to raise warning about slab size.
    
    Returns: 
        list of dictionaries with slabs, list of repeat slabs, list of slabs 
        larger than the max_size 
    """
    # Iterate though provisional slabs to extract the unique slabs
    unique_list, unique_list_of_dicts, repeat, large = ([] for i in range(4))

    for slab in provisional:
        if slab['slab'] not in unique_list:
            unique_list.append(slab['slab'])
            unique_list_of_dicts.append(slab)
            
            # For large slab size warning
            atoms = len(slab['slab'].atomic_numbers)
            if atoms > max_size:
                large.append('{}_{}_{}_{}'.format(slab['hkl'], 
                slab['slab_thickness'], slab['vac_thickness'], slab['slab_index']))

        # For repeat slabs warning
        else:
            repeat.append('{}_{}_{}_{}'.format(slab['hkl'], 
            slab['slab_thickness'], slab['vac_thickness'], slab['slab_index']))

    return unique_list_of_dicts, repeat, large    

def _mp_generate_slabs(struc, hkl, thickness, vacuum, is_symmetric=True, 
center_slab=True, **mp_kwargs): 

    """
    Helper function for multiprocessing in ``slabs_get_single_hkl``, made so 
    that the information on slab and vacuum thickness is noted and carried 
    forward during multiprocessing. ``**mp_kwargs`` can take any 
    ``SlabGenerator`` argument.

    Args: 
        struc (`pymatgen Structure object`): Structure object decorated with 
            oxidation states 
        max_index (`int`): The maximum Miller index to go up to.
        thickness (`int`): Minimum slab thickness 
        vacuum (`int`): Minimum vacuum thickness
        is_symmetric (`bool`, optional): Whether the slabs cleaved should 
            have inversion symmetry. If bulk is non-centrosymmetric, 
            ``is_symmetric`` needs to be ``False`` - the function will return no
            slabs as it looks for inversion symmetry. Take care checking the 
            slabs for mirror plane symmetry before just using them. Defaults to 
            ``True``.  
        center_slab (`bool`, optional): The position of the slab in the 
            simulation cell. 
            
            * ``True``: the slab is centered with equal amounts of 
              vacuum above and below.

            * ``False``: the slab is at the bottom of the simulation cell with
              all of the vacuum on top of it. 

            Defaults to True. 

    Returns
        List of dicts of slabs and relevant metadata
    
    """

    SlabGenerator_kwargs = {'in_unit_planes': False, 'primitive': True, 
    'max_normal_search': None, 'reorient_lattice': True, 'lll_reduce': True}
    SlabGenerator_kwargs.update(
        (k, mp_kwargs[k]) for k in SlabGenerator_kwargs.keys() & mp_kwargs.keys()
        )
    
    get_slabs_kwargs = {'ftol': 0.1, 'tol': 0.1, 'max_broken_bonds': 0, 
    'symmetrize': False, 'repair': False, 'bonds': None}
    get_slabs_kwargs.update(
        (k, mp_kwargs[k]) for k in get_slabs_kwargs.keys() & mp_kwargs.keys()
        )

    slabs = []
    slabgen = SlabGenerator(struc, hkl, thickness, vacuum, 
    center_slab=center_slab, **SlabGenerator_kwargs)
    all_slabs = slabgen.get_slabs(**get_slabs_kwargs)

    h = slabgen._proj_height
    p = round(h/slabgen.parent.lattice.d_hkl(slabgen.miller_index), 8)
    if slabgen.in_unit_planes:
        nlayers_slab = int(math.ceil(slabgen.min_slab_size / p))
    else: 
        nlayers_slab = int(math.ceil(slabgen.min_slab_size / h))

    for i, slab in enumerate(all_slabs):
        if is_symmetric == True:
            if not slab.is_polar() and slab.is_symmetric(): 
                slabs.append({
                'hkl': ''.join(map(str, slab.miller_index)),
                'slab_thickness': thickness,
                'slab_layers': nlayers_slab,
                'vac_thickness': vacuum,
                'slab_index': i,
                'slab': slab})
                
        else: 
            if not slab.is_polar(): 
                slabs.append({
                'hkl': ''.join(map(str, slab.miller_index)),
                'slab_thickness': thickness,
                'slab_layers': nlayers_slab,
                'vac_thickness': vacuum,
                'slab_index': i,
                'slab': slab})
    
    return slabs 

def _get_selective_dynamics(structure, slabs, layers_to_relax=None): 

    # get formula and number of atoms in each primitive unit
    formula = structure.get_primitive_structure().formula 
    primitive_els = ''.join([i for i in formula if not i.isdigit()]).split(' ') 
    primitive_els_num = [int(i) for i in formula if i.isdigit()]
    primitive_atoms = dict(zip(primitive_els, primitive_els_num))
    
    small = []
    for slab in slabs: 
        # Use a temporary slab to make sure the structure is ordered in the same 
        # way as the primitive bulk 
        temp_slab = slab['slab'].get_sorted_structure()
        atoms_per_layer_num = int(len(temp_slab) / slab['slab_layers'])
        # get the multiplier - we can have multiple primitive units in one layer 
        # so we want to make sure the right number of atoms are allowed to relax
        prim_in_layer = atoms_per_layer_num/sum(primitive_els_num)

        # Make sure slab has enough layers to constrain the required layers - 
        # at least one layer should be fixed in the middle, otherwise all layers
        # are allowed to relax and no selective dynamics is applied
        if slab['slab_layers'] <= 2*layers_to_relax:
            small.append('{}_{}_{}_{}'.format(slab['hkl'], 
            slab['slab_thickness'], slab['vac_thickness'], slab['slab_index']))
            slab['slab'] = temp_slab

        else: 
            # get list of lists of atoms in the structure, grouped by the atom 
            grouped_atoms_list = []
            for i in primitive_els: 
                atoms = []
                for site in temp_slab: 
                    if i == site.specie.symbol: 
                        atoms.append([i])
                grouped_atoms_list.append(atoms)

            sd = [] 
            for i, lst in zip(primitive_els, grouped_atoms_list): 
                arr = np.zeros((len(lst), 3))
                allowed_relax = int(primitive_atoms[i]*prim_in_layer*layers_to_relax)
                arr[0:allowed_relax,] = [1,1,1]
                arr[-allowed_relax:len(lst),] = [1,1,1]
                sd.append(arr.tolist())
            sd = [item for sublist in sd for item in sublist]
            
            temp_slab.add_site_property('selective_dynamics', sd)

            slab['slab'] = temp_slab

    return slabs, small  