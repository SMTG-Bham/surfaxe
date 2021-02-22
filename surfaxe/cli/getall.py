# Misc 
from argparse import ArgumentParser
import yaml
import os
import warnings 

# Surfaxe 
from surfaxe.generation import get_slabs_max_index

def _oxstates_to_dict(ox): 
    keys, values = ([] for i in range(2))
    for i in ox.split(','): 
        m = i.split(':')
        values.append(float(m.pop(1)))
        keys.append(str(m.pop(0)))

    ox_states_dict = dict(zip(keys,values))
    return ox_states_dict

def _get_parser(): 
    parser = ArgumentParser(
        description="""Generates all unique slabs with specified maximum 
        Miller index, minimum slab and vacuum thicknesses. It includes all 
        combinations for multiple zero dipole symmetric terminations for the 
        same Miller index"""
    )

    parser.add_argument('-s', '--structure',
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('--hkl', type=int, default=1,
    help='The maximum Miller index (default: 1)')
    parser.add_argument('-t', '--thicknesses', nargs='+', type=int,
    help='The minimum size of the slab in Angstroms.')
    parser.add_argument('-v', '--vacuums', nargs='+', type=int,
    help='The minimum size of the vacuum in Angstroms.')
    parser.add_argument('-r', '--fols', default=False, action='store_true', 
    help=('Makes folders for each termination and slab/vacuum thickness ' 
          'combinations containing POSCARs (default: False)'))
    parser.add_argument('-f', '--files', default=False, action='store_true', 
    help='Makes INCAR, POTCAR and KPOINTS files in each folder (default: False)')
    parser.add_argument('--max-size', default=500, dest='max_size', type=int,
    help=('The maximum number of atoms in the slab specified to raise warning ' 
          'about slab size. Even if the warning is raised, it still outputs ' 
          'the slabs regardless. (default: 500)'))
    parser.add_argument('--no-center-slab', default=True, dest='center_slab',
    action='store_false', help=('The position of the slab in the simulation cell. ' 
    'Centers slab by default'))
    parser.add_argument('--oxstates-list', default=None, dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('--oxstates-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', help=('Add oxidation states to the structure as ' 
    'a dictionary e.g. "Fe:3,O:-2"'))
    parser.add_argument('--no-save', default=True, action='store_false', 
    dest='save_slabs', 
    help='Whether to save the slabs to file (default: True)')
    parser.add_argument('--no-sym', default=True, action='store_false', 
    dest='sym', help=('Whether the slabs cleaved should have inversion symmetry. '
        'By default searches for slabs with inversion symmetry'))
    parser.add_argument('--fmt', default='poscar',  
    help='Format of output files (default: poscar)')
    parser.add_argument('--name', default='POSCAR',  
    help='Name of the surface slab structure file created (default: POSCAR)')
    parser.add_argument('--config-dict', default='PBEsol_config.json', 
    dest='config_dict', 
    help='Specifies the dictionary used for the generation of the input files')
    parser.add_argument('-i', '--incar', default=None,
    help='Overrides the default INCAR parameter settings')
    parser.add_argument('-k', '--kpoints', default=None,
    help='Overrides the default KPOINTS settings.')
    parser.add_argument('-p', '--potcar', default=None,
    help='Overrides the default POTCAR settings')
    parser.add_argument('--yaml', default=False, action='store_true', 
    help='Read optional args from surfaxe_config.yaml file.')

    return parser

def main(): 
    args = _get_parser().parse_args()

    if args.yaml==True: 
        with open('surfaxe_config.yaml', 'r') as y: 
            yaml_args = yaml.load(y)
        args.update(
            (k, yaml_args[k]) for k in args.keys() and yaml_args.keys()
        ) 

    if args.ox_states_dict: 
        ox_states = args.ox_states_dict 
    elif args.ox_states_list: 
        ox_states = map(float, args.ox_states_list.strip('[]').split(','))
    else: 
        ox_states=None 

    get_slabs_max_index(args.structure, args.hkl, args.thicknesses, args.vacuums, 
    make_fols=args.fols, make_input_files=args.files, max_size=args.max_size, 
    center_slab=args.center_slab, ox_states=ox_states, 
    save_slabs=args.save_slabs, is_symmetric=args.sym, fmt=args.fmt, 
    name=args.name, config_dict=args.config_dict, 
    user_incar_settings=args.incar, user_potcar_settings=args.potcar, 
    user_kpoints_settings=args.kpoints)

if __name__ == "__main__":
    main()