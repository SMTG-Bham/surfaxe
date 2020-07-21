
from pymatgen.core.surface import SlabGenerator, generate_all_slabs
from pymatgen import Structure
from pymatgen.io.vasp.sets import DictSet
import os

PBEsol_slab_config = {
  'INCAR': {'ALGO': 'Normal', 'ADDGRID': False, 'LASPH': True, 'EDIFFG': -0.01, 'EDIFF': 1e-06, 'ENCUT': 500, 'ISIF': 2, 'ISMEAR': 0, 'GGA': 'PS', 'LASPH': True, 'LREAL': 'auto', 'LORBIT': 11, 'LCHARG': False, 'LWAVE': False, 'NELM': 150, 'NSW': 0, 'PREC': 'Accurate', 'SIGMA': 0.01, 'ISYM': 2, 'IWAVPR': 1},
  'KPOINTS': {'reciprocal_density': 90},
  'POTCAR': {'Ac': 'Ac', 'Ag': 'Ag', 'Al': 'Al', 'Ar': 'Ar', 'As': 'As', 'Au':  'Au', 'B': 'B',
             'Ba': 'Ba_sv', 'Be': 'Be_sv', 'Bi': 'Bi', 'Br': 'Br', 'C': 'C', 'Ca': 'Ca_sv',
             'Cd': 'Cd', 'Ce': 'Ce', 'Cl': 'Cl', 'Co': 'Co', 'Cr': 'Cr_pv', 'Cs': 'Cs_sv', 'Cu': 'Cu', 'Dy': 'Dy_3', 'Er': 'Er_3', 'Eu': 'Eu', 'F': 'F', 'Fe': 'Fe_pv', 'Ga': 'Ga_d', 'Gd': 'Gd', 'Ge': 'Ge_d', 'H': 'H', 'He': 'He','Hf': 'Hf_pv', 'Hg': 'Hg', 'Ho': 'Ho_3', 'I': 'I', 'In': 'In_d', 'Ir': 'Ir', 'K': 'K_sv', 'Kr': 'Kr', 'La': 'La', 'Li': 'Li_sv', 'Lu': 'Lu_3', 'Mg': 'Mg_pv', 'Mn': 'Mn_pv', 'Mo': 'Mo_pv', 'N': 'N', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'Ne': 'Ne', 'Ni': 'Ni_pv', 'Np': 'Np', 'O': 'O', 'Os': 'Os_pv', 'P': 'P', 'Pa': 'Pa', 'Pb': 'Pb_d', 'Pd': 'Pd', 'Pm': 'Pm_3', 'Pr': 'Pr_3', 'Pt': 'Pt', 'Pu': 'Pu', 'Rb': 'Rb_sv', 'Re': 'Re_pv', 'Rh': 'Rh_pv', 'Ru': 'Ru_pv', 'S': 'S', 'Sb': 'Sb', 'Sc': 'Sc_sv', 'Se': 'Se', 'Si': 'Si', 'Sm': 'Sm_3', 'Sn': 'Sn_d', 'Sr': 'Sr_sv', 'Ta': 'Ta_pv', 'Tb': 'Tb_3', 'Tc': 'Tc_pv', 'Te': 'Te', 'Th': 'Th', 'Ti': 'Ti_sv', 'Tl': 'Tl_d', 'Tm': 'Tm_3', 'U': 'U', 'V': 'V_pv', 'W': 'W', 'Xe': 'Xe', 'Y': 'Y_sv', 'Yb': 'Yb_2', 'Zn': 'Zn', 'Zr': 'Zr_sv'}}

def get_one_hkl_slab(structure, hkl, thicknesses, vacuums, make_fols=False, make_input_files=False, lll_reduce=True, center_slab=True, ox_states=None, max_size=500, potcar_functional='PBE', update_incar=None, update_potcar=None, update_kpoints=None, **kwargs):
    """
    generates all unique slabs for a specified hkl with min slab and vacuum thicknesses
    Args:
        structure (str): filename of structure file
        hkl (tuple): miller index
        thicknesses (list): min size in angstroms of the slab
        vacuums (list): min size in angstroms of vacuum
        ox_states (dict): dict of oxidation states. E.g., {“Li”:1, “Fe”:2, “P”:5, “O”:-2}, defaults to None
        lll_reduce (bool): if true it keeps the slab angles reasonable-ish, default=True
        center_slab (bool): if true slab is centred to the middle of unit cell with equal layers of vacuum above and below, default=True
        update_incar (dict): overrides default INCAR settings; default=None
        update_kpoints (dict or kpoints object): overrides default kpoints settings, if supplied as dict should be as {'reciprocal_density': 100}; default=None
        update_potcar (dict): overrides default POTCAR settings; default=None
        **kwargs allow all SlabGenerator and DictSet args to be passed to the functions
    Returns:
        POSCAR_hkl_slab_vac_index.vasp or hkl/slab_vac_index folders with POSCARs or hkl/slab_vac_index with all input files
    """
    #import bulk relaxed structure, add oxidation states for slab dipole calculations
    struc = Structure.from_file(structure)
    if ox_states is not None:
        struc.add_oxidation_state_by_element(ox_states)
    else:
        struc.add_oxidation_state_by_guess()

    #iterate through vacuums and thicknessses to get all zero dipole symmetric slabs
    provisional = []
    for thickness in thicknesses:
        for vacuum in vacuums:
            slabgen = SlabGenerator(struc, hkl, thickness, vacuum,lll_reduce=lll_reduce, center_slab=center_slab)
            slabs = slabgen.get_slabs()
            for i, slab in enumerate(slabs):
                if (slab.is_polar() == False) and (slab.is_symmetric() == True):
                    provisional.append({'hkl': ''.join(map(str, slab.miller_index)), 'slab_t': thickness, 'vac_t': vacuum, 's_index': i, 'slab': slab})

    #iterate though provisional slabs to extract the unique slabs
    unique_list = []
    unique_list_of_dicts = []
    repeat = []
    large = []

    for slab in provisional:
        if slab['slab'] not in unique_list:
            unique_list.append(slab['slab'])
            unique_list_of_dicts.append(slab)
            atoms = len(slab['slab'].atomic_numbers)
            if atoms > max_size:
                large.append('{}_{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index']))
        else:
            repeat.append('{}_{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index']))

    #warnings for large and repeated slabs
    if repeat:
        print('warning: not all combinations slab/vac thickness were generated because of repeat structures. make sure to manually check the missing slabs.')
        print('the repeat slabs are: ' + ', '.join(map(str, repeat)))

    if large:
        print('warning: some generated slabs exceed the max size specified')
        print('slabs that exceed the max size are: ' + ', '.join(map(str, large)))

    #makes folders hkl/slab_vac_index with POSCAR, POTCAR, INCAR and KPOINTS
    if make_fols is True:
        for root, fols, files in os.walk(os.getcwd()):
            for slab in unique_list_of_dicts:
                        os.mkdir(os.path.join(root, r'{}/{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index'])))
                        if make_input_files is True:
                            vis = DictSet(structure=slab['slab'], config_dict=PBEsol_slab_config, potcar_functional=potcar_functional, user_incar_settings=update_incar, user_potcar_settings=update_potcar, user_kpoints_settings=update_kpoints)
                            vis.write_input(os.path.join(root,r'{}/{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index'])))
                        else:
                            slab['slab'].to(fmt='poscar',filename=r'{}/{}_{}_{}/POSCAR'.format(slab['hkl'],slab['slab_t'], slab['vac_t'], slab['s_index']))

    #omits folders, makes POSCAR_hkl_slab_vac_index files in the root directory
    else:
        for slab in unique_list_of_dicts:
            slab['slab'].to(fmt='poscar',filename='POSCAR_{}_{}_{}_{}.vasp'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index']))



def get_all_unique_slabs(structure, max_index, thicknesses, vacuums, make_fols=False, make_input_files=False, max_size=500, ox_states=None, lll_reduce=True, center_slab=True, potcar_functional='PBE', update_incar=None, update_potcar=None, update_kpoints=None, **kwargs):
    """
    generates all unique slabs with specified max hkl, min slab and vacuum thicknesses; including all combinations for multiple zero-dipole symmetric terminations for the same hkl
    Args:
        structure: filename of structure file; oxidation states added by guess
        max_index (int): maximum Miller index to be considered
        thicknesses (list): min size in of the slab in angstroms
        vacuums (list): min size in of vacuum in angstroms
        make_fols (bool): makes folders containing POSCARs, default=False
        make_input_files (bool): makes INCAR, POTCAR and KPOINTS files in each of the folders, default=False
        max_size (int): allows a warning based on the number of atoms, default=500
        ox_states (dict): dict of oxidation states. E.g., {“Li”:1, “Fe”:2, “P”:5, “O”:-2}, default=None
        lll_reduce (bool): if true it keeps the slab angles reasonable-ish, default=True
        center_slab (bool): slab is centred to the middle of unit cell with equal layers of vacuum above and below, default=True
        update_incar (dict): overrides default INCAR settings; default=None
        update_kpoints (dict or kpoints object): overrides default kpoints settings, if supplied as dict should be as {'reciprocal_density': 100}; default=None
        update_potcar (dict): overrides default POTCAR settings; default=None
        **kwargs allow all SlabGenerator and DictSet args to be passed to the functions
    Returns:
        POSCAR_hkl_slab_vac_index.vasp or hkl/slab_vac_index folders with POSCARs or hkl/slab_vac_index with all input files

    """
    #import bulk relaxed structure, add oxidation states for slab dipole calculations
    struc = Structure.from_file(structure)
    if ox_states is not None:
        struc.add_oxidation_state_by_element(ox_states)
    else:
        struc.add_oxidation_state_by_guess()

    #iterate through vacuums and thicknessses to get all zero dipole symmetric slabs
    provisional = []
    for vacuum in vacuums:
        for thickness in thicknesses:
            all_slabs = generate_all_slabs(struc, max_index=max_index, min_slab_size=thickness, min_vacuum_size=vacuum, lll_reduce=lll_reduce, center_slab=center_slab)
            for i, slab in enumerate(all_slabs):
                if (slab.is_polar() == False) and (slab.is_symmetric() == True):
                    provisional.append({'hkl': ''.join(map(str, slab.miller_index)), 'slab_t': thickness, 'vac_t': vacuum, 's_index': i, 'slab': slab})

    #iterate though provisional slabs to extract the unique slabs
    unique_list = []
    unique_list_of_dicts = []
    repeat = []
    large = []

    for slab in provisional:
        if slab['slab'] not in unique_list:
            unique_list.append(slab['slab'])
            unique_list_of_dicts.append(slab)
            atoms = len(slab['slab'].atomic_numbers)
            if atoms > max_size:
                large.append('{}_{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index']))

        else:
            repeat.append('{}_{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index']))

    #warnings for large and repeated slabs
    if repeat:
        print('warning: not all combinations of hkl or slab/vac thickness were generated because of repeat structures.')
        print('the repeat slabs are: ' + ', '.join(map(str, repeat)))

    if large:
        print('warning: some generated slabs exceed the max size specified')
        print('slabs that exceed the max size are: ' + ', '.join(map(str, large)))

    #makes folders hkl/slab_vac_index with POSCAR, POTCAR, INCAR and KPOINTS
    if make_fols is True:
        for root, fols, files in os.walk(os.getcwd()):
            for slab in unique_list_of_dicts:
                os.mkdir(os.path.join(root, r'{}/{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index'])))
                if make_input_files is True:
                    vis = DictSet(structure=slab['slab'],config_dict=PBEsol_slab_config, potcar_functional=potcar_functional, user_incar_settings=update_incar, user_potcar_settings=update_potcar, user_kpoints_settings=update_kpoints)
                    vis.write_input(os.path.join(root, r'{}/{}_{}_{}'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index'])))
                else:
                    slab['slab'].to(fmt='poscar',filename=r'{}/{}_{}_{}/POSCAR'.format(slab['hkl'],slab['slab_t'], slab['vac_t'], slab['s_index']))

    #omits folders, makes POSCAR_hkl_slab_vac_index files in the root folder
    else:
        for slab in unique_list_of_dicts:
            slab['slab'].to(fmt='poscar',filename='POSCAR_{}_{}_{}_{}.vasp'.format(slab['hkl'], slab['slab_t'], slab['vac_t'], slab['s_index']))
