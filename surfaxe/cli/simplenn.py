# Misc 
from argparse import ArgumentParser
import json
import os
import warnings 
from pymatgen.analysis.local_env import CrystalNN

# Surfaxe 
from surfaxe.analysis import simple_nn

def _get_parser(): 
    parser = ArgumentParser(
        description="""Finds the nearest neighbours for simple structures. 
        Before using on slabs make sure the nn_method works with the bulk 
        structure."""
    )
    
    parser.add_argument('-s', '--start', required=True, type=str, 
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('-a', '--atoms', required=True, 
    help='List of elements in the structure in any order')
    parser.add_argument('-e', '--end', type=str, default=None,
    help=('Filename of structure file in any format supported by pymatgen. ' 
          'Use if comparing initial and final structures.'))
    parser.add_argument('--oxstates-list', default=None, dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('--oxstates-dict', default=None, type=dict,
    dest='ox_states_dict', 
    help='Add oxidation states to the structure as a dictionary.')
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='nn_data.csv', type=str,
    dest='csv_fname', help='Filename of the csv file (default: nn_data.csv)')
    
    return parser 

def main(): 
    args = _get_parser().parse_args()
    elements = map(str, args.atoms.strip('[]').split(','))

    if args.ox_states_dict: 
        ox_states = args.ox_states_dict 
    elif args.ox_states_list: 
        ox_states = map(float, args.ox_states_list.strip('[]').split(','))
    else: 
        ox_states=None
    
    if args.save_csv is True: 
        simple_nn(args.start, elements, end=args.end, ox_states=ox_states, 
        nn_method=CrystalNN(), save_csv=args.save_csv, csv_fname=args.csv_fname)
    else: 
        nn = simple_nn(args.start, elements, end=args.end, ox_states=ox_states, 
        nn_method=CrystalNN(), save_csv=args.save_csv, csv_fname=args.csv_fname)
        print(nn)

if __name__ == "__main__":
    main()