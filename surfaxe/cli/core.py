# Misc 
from argparse import ArgumentParser
import os
import yaml
from pymatgen.analysis.local_env import CrystalNN

# Surfaxe 
from surfaxe.vasp_data import core_energy

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
        description="""Parses the structure and OUTCAR files for the core level 
        energy. Check the validity of nearest neighbour method on the bulk 
        structure before using it on slabs."""
    )

    parser.add_argument('-a', '--atom', dest='core_atom', 
    help='The symbol of atom the core state energy level should be parsed from')
    parser.add_argument('-b', '--bulknn',  nargs='+', type=str,
    help=('The symbols of the nearest neighbours of the core atom ' 
    'e.g. "Ti Ti Y Y"'))
    parser.add_argument('-o', '--orbital', default='1s', 
    help='The orbital of core state (default: 1s)')
    parser.add_argument('--oxstates-list', default=None, dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('--oxstates-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', help=('Add oxidation states to the structure as ' 
    'a dictionary e.g. "Fe:3,O:-2"'))
    parser.add_argument('-s', '--structure', default='POSCAR',
    help=('Filename of structure file in any format supported by pymatgen '
          '(default: POSCAR'))
    parser.add_argument('--outcar', default='OUTCAR', 
    help='Path to OUTCAR file (default: OUTCAR)')
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

    core = core_energy(args.core_atom, args.bulk_nn, orbital=args.orbital, 
    ox_states=ox_states, nn_method=CrystalNN(), structure=args.structure)
    print(core)
    
if __name__ == "__main__":
    main()