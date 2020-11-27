# Misc 
from argparse import ArgumentParser
import json
import os
import warnings 

# Surfaxe 
from surfaxe.analysis import bond_analysis

def _get_parser(): 
    parser = ArgumentParser(
        description="""Parses the structure looking for bonds between atoms. 
        Check the validity of the nearest neighbour method on the bulk structure 
        before using it on slabs."""
    )

    parser.add_argument('-s', '--structure', required=True, type=str,
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('-b', '--bonds', required=True, type=list, 
    help='List of bonds to compare in any order')
    # not sure how to pass class as the default here? 
    parser.add_argument('-n', '--nnmethod', default='CrystalNN()', type=str, 
    dest='nn_method', 
    help='The pymatgen local_env nearest neighbour method (default: CrystalNN()')
    parser.add_argument('--oxstates', default=None, dest='ox_states',
    help='Add oxidation states to the structure.')
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='bond_analysis.csv', type=str,
    dest='csv_fname', help='Filename of the csv file (default: bond_analysis.csv)')
    parser.add_argument('--no-plot', default=True, action='store_false', 
    dest='save_plt', help='Turns off plotting the bond lengths')
    parser.add_argument('--plt-fname', default='bond_analysis.png', type=str,
    dest='plt_fname', help='Filename of the plot')
    parser.add_argument('--dpi', default=300, type=int, help='Dots per inch')

    return parser 

def main(): 
    args = _get_parser().parse_args()

    # warnings? 
    
    bond_analysis(args.structure, args.bonds, nn_method=args.nn_method, 
    ox_states=args.ox_states, save_csv=args.save_csv, csv_fname=args.csv_fname, 
    save_plt=args.save_plt, plt_fname=args.plt_fname, dpi=args.dpi)

if __name__ == "__main__":
    main()