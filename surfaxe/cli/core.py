# Misc 
from argparse import ArgumentParser
import os
import yaml
from pymatgen.analysis.local_env import CrystalNN

# Surfaxe 
from surfaxe.data import core_energy

def _get_parser(): 
    parser = ArgumentParser(
        description="""Parses the structure and OUTCAR files for the core level 
        energy. Check the validity of nearest neighbour method on the bulk 
        structure before using it on slabs."""
    )

    parser.add_argument('-p', '--path', default=None, help=('The path to the ' 
    'structure and OUTCAR files'))
    parser.add_argument('-a', '--atom', dest='core_atom', 
    help='The symbol of atom the core state energy level should be parsed from')
    parser.add_argument('-b', '--bulknn',  
    help='The symbols of the nearest neighbours of the core atom e.g. "Ti,Ti,Y,Y"')
    parser.add_argument('-o', '--orbital', default='1s', 
    help='The orbital of core state (default: 1s)')
    parser.add_argument('--oxstates-list', default=None, dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('-s', '--structure', default='POSCAR', 
    help=('Filename of structure file in any format supported by pymatgen '
    'default: POSCAR'))
    parser.add_argument('--yaml', default=False, action='store_true', 
    help='Read optional args from surfaxe_config.yaml file.')

    return parser

def main(): 
    args = _get_parser().parse_args()
    args.ox_states_dict = None
    bulk_nn = map(str, args.bulknn.strip('[]').split(','))

    path = os.getcwd()
    if args.path is not None: 
        path = args.path

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

    core = core_energy(path, args.core_atom, bulk_nn, orbital=args.orbital, 
    ox_states=ox_states, nn_method=CrystalNN(), structure=args.structure)
    print(core)
    
if __name__ == "__main__":
    main()