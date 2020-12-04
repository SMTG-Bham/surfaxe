# Misc 
from argparse import ArgumentParser
import json
import os
import warnings 

# Surfaxe 
from surfaxe.analysis import complex_nn

def _get_parser(): 
    parser = ArgumentParser(
        description="""Finds the nearest neighbours for more complex structures. 
        Uses CutOffDictNN() class as the nearest neighbour method. Check 
        validity on bulk structure before applying to surface slabs."""
    )
    
    parser.add_argument('-s', '--start', required=True, type=str, 
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('-a', '--atoms', required=True, type=list, dest='elements',
    help='List of elements in the structure in any order')
    parser.add_argument('-c', '--cutoffdict', required=True, type=dict, 
    dest='cut_off_dict', help='Dictionary of bond lengths')
    parser.add_argument('-e', '--end', type=str, default=None,
    help=('Filename of structure file in any format supported by pymatgen. ' 
          'Use if comparing initial and final structures.'))
    parser.add_argument('--oxstates-list', default=None, type=list,
    dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('--oxstates-dict', default=None, type=dict,
    dest='ox_states_dict', 
    help='Add oxidation states to the structure as a dictionary.')
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='nn_data.csv', type=str,
    dest='csv_fname', help='Filename of the csv file (default: nn_data.csv)')
    
    return parser 

def _ox_state_formatter(): 
    pass

def main(): 
    args = _get_parser().parse_args()

    if args.ox_states_dict: 
        ox_states = args.ox_states_dict 
    elif args.ox_states_list: 
        ox_states = args.ox_states_list
    else: 
        ox_states=None 

    complex_nn(args.start, args.elements, args.cut_off_dict, end=args.end, 
    ox_states=ox_states, save_csv=args.save_csv, csv_fname=args.csv_fname)

if __name__ == "__main__":
    main()