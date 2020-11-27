# Misc 
from argparse import ArgumentParser
import json
import os
import warnings 

# Surfaxe 
from surfaxe.analysis import simple_nn

def _get_parser(): 
    parser = ArgumentParser(
        description="""Finds the nearest neighbours for simple structures. 
        Before using on slabs make sure the nn_method works with the bulk 
        structure."""
    )
    
    parser.add_argument('-s' '--start', required=True, type=str, 
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('-a', '--atoms', required=True, type=list, dest='elements',
    help='List of elements in the structure in any order')
    parser.add_argument('-e', '--end', type=str, default=None,
    help=('Filename of structure file in any format supported by pymatgen. ' 
          'Use if comparing initial and final structures.'))
    # not sure what to do if theres more than one supported type 
    parser.add_argument('--oxstates', default=None, dest='ox_states',
    help='Add oxidation states to the structure.')
    # not sure how to pass class as the default here? 
    parser.add_argument('-n', '--nnmethod', default='CrystalNN()', type=str, 
    dest='nn_method'
    help='The pymatgen local_env nearest neighbour method (default: CrystalNN()')
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='nn_data.csv', type=str,
    dest='csv_fname', help='Filename of the csv file (default: nn_data.csv)')
    
    return parser 

def main(): 
    args = _get_parser().parse_args()

    # warnings? 

    simple_nn(args.start, args.elements, end=args.end, ox_states=args.ox_states, 
    nn_method=args.nn_method, save_csv=args.save_csv, csv_fname=args.csv_fname)

if __name__ == "__main__":
    main()