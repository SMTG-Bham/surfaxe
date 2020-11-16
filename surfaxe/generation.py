
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen import Structure
from pymatgen.io.vasp.sets import DictSet
import warnings
import os

# surfaxe
from surfaxe.io import slabs_to_file, _custom_formatwarning

def get_one_hkl_slabs(structure, hkl, thicknesses, vacuums, make_fols=False, 
make_input_files=False, max_size=500, bonds=None, center_slab=True, 
ox_states=None, save_slabs=True, is_symmetric=True, 
config_dict='PBEsol_config.json', user_incar_settings=None, 
user_kpoints_settings=None, user_potcar_settings=None, **kwargs):
    """
    Generates all unique slabs for a specified Miller index with minimum slab
    and vacuum thicknesses. 
    
    Note that using this method of slab generation will result in different slab 
    index numbers as in the `get_all_slabs` - the slabs identified are the same, 
    the index varies based on the position in the list of generated slabs. 
    The function returns None by default and generates either: 

    (i) POSCAR_hkl_slab_vac_index.vasp (default) 
    (ii) hkl/slab_vac_index folders with POSCARs
    (iii) hkl/slab_vac_index with all VASP input files 
    
    Or if `save_slabs=False` a list of dicts of all unique slabs is returned. 
    
    Args:
        structure (`str`): Filename of structure file in any format supported by 
            pymatgen. 
        hkl (`tuple`): Miller index. Defaults to ``None``.
        thicknesses (`list`): The minimum size of the slab in Angstroms. 
            Defaults to ``None``. 
        vacuums (`list`): The minimum size of the vacuum in Angstroms. Defaults 
            to ``None``. 
        make_fols (`bool`, optional): Makes folders for each termination 
            and slab/vacuum thickness combinations containing POSCARs. 
            
            * ``True``: A Miller index folder is created, in which folders 
              named slab_vac_index are created to which the relevant POSCARs 
              are saved. 
                    
                    E.g. for a (0,0,1) slab of index 1 with a slab thickness of 
                    20 Å and vacuum thickness of 30 Å the folder structure would 
                    be: ``001/20_30_1/POSCAR``  

            * ``False``: The indexed POSCARs are put in a folder named after 
              the bulk formula. 
              
                    E.g. for a (0,0,1) MgO slab of index 1 with a slab thickness 
                    of 20 Å and vacuum thickness of 30 Å the folder structure 
                    would be: ``MgO/POSCAR_001_20_30_1``

            Defaults to ``False``.  
        make_input_files (`bool`, optional): Makes INCAR, POTCAR and 
            KPOINTS files in each folder. If ``make_input_files`` is ``True`` 
            but ``make_files`` or ``save_slabs`` is ``False``, files will be 
            saved to folders regardless. Defaults to ``False``. 
        max_size (`int`, optional): The maximum number of atoms in the slab 
            specified to raise warning about slab size. Even if the warning is 
            raised, it still outputs the slabs regardless. Defaults to ``500``. 
        bonds ({(specie1, specie2): max_bond_dist}: bonds are specified as  
            {string tuple: float} of specie1, specie2 and the max bonding 
            distance. For example, PO4 groups may be defined as {(“P”, “O”): 3}.
        center_slab (`bool`, optional): The position of the slab in the 
            simulation cell. 
            
            * ``True``: the slab is centered with equal amounts of 
              vacuum above and below.

            * ``False``: the slab is at the bottom of the simulation cell with
              all of the vacuum on top of it. 

            Defaults to True. 

        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess. 
              
            Defaults to ``None``. 

        save_slabs (`bool`, optional): Whether to save the slabs to file. 
            Defaults to ``True``.
        is_symmetric (`bool`, optional): Whether the slabs cleaved should 
            have inversion symmetry. If bulk is non-centrosymmetric, 
            ``is_symmetric`` needs to be ``False`` - the function will return no
            slabs as it looks for inversion symmetry. Take care checking the 
            slabs for mirror plane symmetry before just using them. Defaults to 
            ``True``. 
        config_dict (`dict` or `str`, optional): Specifies the dictionary used 
            for the generation of the input files. Defaults to ``PBEsol_config.json`` 
        user_incar_settings (`dict`, optional): Overrides the default INCAR 
            parameter settings. Defaults to ``None``.
        user_kpoints_settings (`dict` or Kpoints object, optional): 
            Overrides the default kpoints settings. If it is supplied  
            as `dict`, it should be as ``{'reciprocal_density': 100}``. Defaults 
            to ``None``.
        user_potcar_settings (`dict`, optional): Overrides the default POTCAR 
            settings. Defaults to ``None``.

    Returns:
        Surface slabs 
    """

    # Set up additional arguments for slab generation and saving slabs
    SlabGenerator_kwargs = {'in_unit_planes': False, 'primitive': True, 
    'max_normal_search': None, 'reorient_lattice': True, 'lll_reduce': True}
    SlabGenerator_kwargs.update(
        (k, kwargs[k]) for k in SlabGenerator_kwargs.keys() & kwargs.keys()
        )
    
    get_slabs_kwargs = {'ftol': 0.1, 'tol': 0.1, 'max_broken_bonds': 0, 
    'symmetrize': False, 'repair': False}
    get_slabs_kwargs.update(
        (k, kwargs[k]) for k in get_slabs_kwargs.keys() & kwargs.keys()
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

    # Import bulk relaxed structure, add oxidation states for slab dipole
    # calculations
    struc = Structure.from_file(structure)
    struc = oxidation_states(struc, ox_states=ox_states)

    # Iterate through vacuums and thicknessses 
    provisional = []
    for vacuum in vacuums:
        for thickness in thicknesses:
            slabgen = SlabGenerator(struc, hkl, thickness, vacuum,
                                    center_slab=center_slab,  
                                    **SlabGenerator_kwargs) 
                                    
            slabs = slabgen.get_slabs(bonds, **get_slabs_kwargs)
            for i, slab in enumerate(slabs):
                # Get all the zero-dipole slabs with inversion symmetry
                if is_symmetric: 
                    if slab.is_symmetric() and not slab.is_polar():
                        provisional.append({
                            'hkl': ''.join(map(str, slab.miller_index)),
                            'slab_t': thickness,
                            'vac_t': vacuum,
                            's_index': i,
                            'slab': slab})
                
                # Get all the zero-dipole slabs wihtout inversion symmetry
                else: 
                    if not slab.is_polar():
                        provisional.append({
                            'hkl': ''.join(map(str, slab.miller_index)),
                            'slab_t': thickness,
                            'vac_t': vacuum,
                            's_index': i,
                            'slab': slab})
                  
    ### DWD: This is block is repeated almost exactly in both the functions
    ### so could we farm it out to another function e.g. filter_slabs() 
    ### so that it is only written once and makes the two main functions
    ### a bit shorter? 
    # Iterate though provisional slabs to extract the unique slabs
    unique_list, unique_list_of_dicts, repeat, large = ([] for i in range(4))

    for slab in provisional:
        if slab['slab'] not in unique_list:
            unique_list.append(slab['slab'])
            unique_list_of_dicts.append(slab)
            # For large slab size warning
            atoms = len(slab['slab'].atomic_numbers)
            if atoms > max_size:
                large.append('{}_{}_{}_{}'.format(slab['hkl'], slab['slab_t'],
                             slab['vac_t'], slab['s_index']))

        # For repeat slabs warning
        else:
            repeat.append('{}_{}_{}_{}'.format(slab['hkl'], slab['slab_t'],
                          slab['vac_t'], slab['s_index']))

    # Warnings for large and repeated slabs
    if repeat:
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Not all combinations of hkl or slab/vac thicknesses '
        'were generated because of repeat structures. '
        'The repeat slabs are: ' + ', '.join(map(str, repeat)))

    if large:
        warnings.formatwarning = _custom_formatwarning
        warnings.warn('Some generated slabs exceed the max size specified.'
        ' Slabs that exceed the max size are: ' + ', '.join(map(str, large)))

    # Save the slabs to file or return the list of dicts 
    if save_slabs: 
        slabs_to_file(list_of_slabs=unique_list_of_dicts, structure=structure, 
        make_fols=make_fols, make_input_files=make_input_files, 
        config_dict=config_dict, **save_slabs_kwargs)
    
    else: 
        return unique_list_of_dicts

def get_all_slabs(structure, max_index, thicknesses, vacuums, make_fols=False, 
make_input_files=False, max_size=500, bonds=None, center_slab=True, 
ox_states=None, save_slabs=True, is_symmetric=True, config_dict=None, 
user_incar_settings=None, user_potcar_settings=None, user_kpoints_settings=None, 
**kwargs):
    """
    Generates all unique slabs with specified maximum Miller index, minimum slab
    and vacuum thicknesses. It includes all combinations for multiple zero
    dipole symmetric terminations for the same Miller index. 
    
    Note that using this method of slab generation will results in different 
    slab index values as in the `get_one_hkl_slabs` - the slabs identified are 
    the same, the index varies based on the position in the list of generated 
    slabs.
    The function returns None by default and generates either: 

    (i) POSCAR_hkl_slab_vac_index.vasp (default) 
    (ii) hkl/slab_vac_index folders with POSCARs
    (iii) hkl/slab_vac_index with all VASP input files 
    
    Or if `save_slabs=False` a list of dicts of all unique slabs is returned. 

    Args:
        structure (`str`, required): Filename of structure file in any 
            format supported by pymatgen. 
        max_index (`int`, required): The maximum Miller index to go up to.
        thicknesses (`list`, required): The minimum size of the slab in 
            Angstroms. 
        vacuums (`list`, required): The minimum size of the vacuum in 
            Angstroms. 
        make_fols (`bool`, optional): Makes folders for each termination 
            and slab/vacuum thickness combinations containing POSCARs. 
            
            * ``True``: A Miller index folder is created, in which folders 
              named slab_vac_index are created to which the relevant POSCARs 
              are saved. 
                    
                    E.g. for a (0,0,1) slab of index 1 with a slab thickness of 
                    20 Å and vacuum thickness of 30 Å the folder structure would 
                    be: ``001/20_30_1/POSCAR``  

            * ``False``: The indexed POSCARs are put in a folder named after 
              the bulk formula. 
              
                    E.g. for a (0,0,1) MgO slab of index 1 with a slab thickness 
                    of 20 Å and vacuum thickness of 30 Å the folder structure 
                    would be: ``MgO/POSCAR_001_20_30_1``

            Defaults to ``False``.  
        make_input_files (`bool`, optional): Makes INCAR, POTCAR and 
            KPOINTS files in each folder. If ``make_input_files`` is ``True`` 
            but ``make_files`` or ``save_slabs`` is ``False``, files will be 
            saved to folders regardless. Defaults to ``False``. 
        max_size (`int`, optional): The maximum number of atoms in the slab 
            specified to raise warning about slab size. Even if the warning is 
            raised, it still outputs the slabs regardless. Defaults to ``500``. 
        bonds ({(specie1, specie2): max_bond_dist}: bonds are specified as  
            {string tuple: float} of specie1, specie2 and the max bonding 
            distance. For example, PO4 groups may be defined as {(“P”, “O”): 3}.
        center_slab (`bool`, optional): The position of the slab in the 
            simulation cell. 
            
            * ``True``: the slab is centered with equal amounts of 
              vacuum above and below.

            * ``False``: the slab is at the bottom of the simulation cell with
              all of the vacuum on top of it. 

            Defaults to True. 

        ox_states (``None``, `list` or  `dict`, optional): Add oxidation states 
            to the structure. Different types of oxidation states specified will 
            result in different pymatgen functions used. The options are: 
            
            * if supplied as ``list``: The oxidation states are added by site 
                    
                    e.g. ``[3, 2, 2, 1, -2, -2, -2, -2]``
            
            * if supplied as ``dict``: The oxidation states are added by element
                    
                    e.g. ``{'Fe': 3, 'O':-2}``
            
            * if ``None``: The oxidation states are added by guess. 
              
            Defaults to ``None``. 

        save_slabs (`bool`, optional): Whether to save the slabs to file. 
            Defaults to True.
        is_symmetric (`bool`, optional): Whether the slabs cleaved should 
            have inversion symmetry. If bulk is non-centrosymmetric, 
            ``is_symmetric`` needs to be ``False`` - the function will return no
            slabs as it looks for inversion symmetry. Take care checking the 
            slabs for mirror plane symmetry before just using them. Defaults to 
            ``True``. 
        config_dict (`dict` or `str`, optional): Specifies the dictionary used 
            for the generation of the input files. Defaults to ``None`` which 
            loads the ``PBEsol_config.json`` file. 
        user_incar_settings (`dict`, optional): Overrides the default INCAR 
            parameter settings. Defaults to ``None``.
        user_kpoints_settings (`dict` or Kpoints object, optional): 
            Overrides the default kpoints settings. If it is supplied  
            as `dict`, it should be as ``{'reciprocal_density': 100}``. Defaults 
            to ``None``.
        user_potcar_settings (`dict`, optional): Overrides the default POTCAR 
            settings. Defaults to ``None``.

    Returns:
        None (default) 
        or unique_slabs (list of dicts) 

    """
   
    # Set up additional arguments for slab generation and saving slabs
    all_slabs_kwargs = {'in_unit_planes': False, 'primitive': True, 
    'max_normal_search': None, 'lll_reduce': True, 'ftol': 0.1, 'tol': 0.1, 
    'max_broken_bonds': 0, 'symmetrize': False, 'repair': False}
    all_slabs_kwargs.update(
        (k, kwargs[k]) for k in all_slabs_kwargs.keys() & kwargs.keys()
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

    # Import bulk relaxed structure, add oxidation states for slab dipole
    # calculations
    struc = Structure.from_file(structure)
    struc = oxidation_states(struc, ox_states=ox_states)
   
    # Iterate through vacuums and thicknessses
    provisional = []
    for vacuum in vacuums:
        for thickness in thicknesses:
            all_slabs = generate_all_slabs(struc, max_index, thickness, vacuum,
                                        center_slab=center_slab, bonds=bonds, 
                                        **all_slabs_kwargs)
            
            for i, slab in enumerate(all_slabs):
                # Get all the zero-dipole slabs with inversion symmetry
                if is_symmetric: 
                    if slab.is_symmetric() and not slab.is_polar():
                        provisional.append({
                            'hkl': ''.join(map(str, slab.miller_index)),
                            'slab_t': thickness,
                            'vac_t': vacuum,
                            's_index': i,
                            'slab': slab})
                
                # Get all the zero-dipole slabs wihtout inversion symmetry 
                else: 
                    if not slab.is_polar():
                        provisional.append({
                            'hkl': ''.join(map(str, slab.miller_index)),
                            'slab_t': thickness,
                            'vac_t': vacuum,
                            's_index': i,
                            'slab': slab
                        })

    # Iterate though provisional slabs to extract the unique slabs
    unique_list, unique_list_of_dicts, repeat, large = ([] for i in range(4))

    for slab in provisional:
        if slab['slab'] not in unique_list:
            unique_list.append(slab['slab'])
            unique_list_of_dicts.append(slab)
            # For large slab size warning
            atoms = len(slab['slab'].atomic_numbers)
            if atoms > max_size:
                large.append('{}_{}_{}_{}'.format(
                    slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index'])
                )

        # For repeat slabs warning
        else:
            repeat.append('{}_{}_{}_{}'.format(
                slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index'])
            )

    # Warnings for large and repeated slabs
        if repeat:
            warnings.formatwarning = _custom_formatwarning
            warnings.warn('Not all combinations of hkl or slab/vac thicknesses '
            'were generated because of repeat structures. '
            'The repeat slabs are: ' + ', '.join(map(str, repeat)))

        if large:
            warnings.formatwarning = _custom_formatwarning
            warnings.warn('Some generated slabs exceed the max size specified.'
            ' Slabs that exceed the max size are: ' + ', '.join(map(str, large)))

    # Save the slabs to file or return the list of dicts 
    if save_slabs: 
        slabs_to_file(list_of_slabs=unique_list_of_dicts, structure=structure, 
        make_fols=make_fols, make_input_files=make_input_files, 
        config_dict=config_dict, **save_slabs_kwargs)
    
    else: 
        return unique_list_of_dicts 

def oxidation_states(structure, ox_states=None):
    ''' 
    Adds oxidation states to the structure object 
    
    Args: 
        structure (`str`, required): Filename of structure file in any 
            format supported by pymatgen. 
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
    if type(ox_states) is dict:
        structure.add_oxidation_state_by_element(ox_states)
    elif type(ox_states) is list:
        structure.add_oxidation_state_by_site(ox_states)
    else:
        structure.add_oxidation_state_by_guess(max_sites=-1)

    return structure