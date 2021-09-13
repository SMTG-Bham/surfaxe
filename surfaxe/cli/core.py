# Misc 
from argparse import ArgumentParser
from ruamel.yaml.main import YAML
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
        structure before using it on slabs. Uses CrystalNN by default."""
    )

    parser.add_argument('-a', '--atom', dest='core_atom', 
    help='The symbol of atom the core state energy level should be parsed from')
    parser.add_argument('--nn',  nargs='+', type=str,
    help=('The symbols of the nearest neighbours of the core atom ' 
    'e.g. "Ti Ti Y Y"'))
    parser.add_argument('-o', '--orbital', default='1s', 
    help='The orbital of core state (default: 1s)')
    parser.add_argument('--oxi-list', default=None, dest='ox_states_list', 
    nargs='+', type=float, 
    help='Add oxidation states to the structure as a list e.g. 3 3 -2 -2 -2')
    parser.add_argument('--oxi-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', help=('Add oxidation states to the structure as ' 
    'a dictionary e.g. Fe:3,O:-2'))
    parser.add_argument('-s', '--structure', default='POSCAR',
    help=('Filename of structure file in any format supported by pymatgen '
          '(default: ./POSCAR'))
    parser.add_argument('--outcar', default='OUTCAR', 
    help='Path to OUTCAR file (default: ./OUTCAR)')
    parser.add_argument('--yaml', default=None, type=str, 
    help=('Read all args from a yaml config file. Completely overrides any '
    'other flags set '))

    return parser

def main(): 
    args = _get_parser().parse_args()
    
    if args.yaml is not None: 
        with open(args.yaml, 'r') as y:  
            yaml = YAML(typ='safe', pure=True)
            yaml_args = yaml.load(y)
        
        core = core_energy(**yaml_args)
        print(core)

    else: 
        if args.ox_states_dict: 
            ox_states = args.ox_states_dict 
        elif args.ox_states_list: 
            ox_states = args.ox_states_list
        else: 
            ox_states=None

        core = core_energy(args.core_atom, args.nn, orbital=args.orbital, 
        ox_states=ox_states, nn_method=CrystalNN(), structure=args.structure)
        print(core)
    
if __name__ == "__main__":
    main()