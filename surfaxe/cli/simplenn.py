# Misc 
from argparse import ArgumentParser
import yaml
import os
import warnings 
from pymatgen.analysis.local_env import CrystalNN

# Surfaxe 
from surfaxe.analysis import simple_nn

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
        description="""Finds the nearest neighbours for simple structures. 
        Before using on slabs make sure the nn_method works with the bulk 
        structure."""
    )
    
    parser.add_argument('-s', '--start', default='POSCAR',
    help=('Filename of structure file in any format supported by pymatgen '
          '(default: POSCAR'))
    parser.add_argument('-a', '--atoms', default=None, nargs='+', type=str,
    help='List of elements in the structure in any order e.g. Y Ti O S')
    parser.add_argument('-e', '--end', default=None,
    help=('Filename of structure file in any format supported by pymatgen. ' 
          'Use if comparing initial and final structures.'))
    parser.add_argument('--oxi-list', default=None, dest='ox_states_list', 
    nargs='+', type=float, 
    help='Add oxidation states to the structure as a list e.g. 3 3 -2 -2 -2')
    parser.add_argument('--oxi-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', help=('Add oxidation states to the structure as ' 
    'a dictionary e.g. Fe:3,O:-2'))
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Prints data to terminal' )
    parser.add_argument('--csv-fname', default='nn_data.csv',
    dest='csv_fname', help='Filename of the csv file (default: nn_data.csv)')
    parser.add_argument('--yaml', default=False, action='store_true', 
    help=('Read all args from surfaxe_config.yaml file. Completely overrides any '
    'other flags set '))
    
    return parser 

def main(): 
    args = _get_parser().parse_args()

    if args.yaml==True: 
        with open('surfaxe_config.yaml', 'r') as y: 
            yaml_args = yaml.safe_load(y)
        
        nn = simple_nn(**yaml_args)
        if ('save_csv', False) in yaml_args.items(): 
            print(nn)

    else: 
        if args.ox_states_dict: 
            ox_states = args.ox_states_dict 
        elif args.ox_states_list: 
            ox_states = args.ox_states_list
        else: 
            ox_states=None

        nn = simple_nn(args.start, args.atoms, end=args.end, ox_states=ox_states, 
        nn_method=CrystalNN(), save_csv=args.save_csv, csv_fname=args.csv_fname)
    
        if args.save_csv==True: 
            print(nn)

if __name__ == "__main__":
    main()