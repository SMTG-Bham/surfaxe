# Misc 
from argparse import ArgumentParser
import json
import os
import warnings 

# Surfaxe 
from surfaxe.convergence import parse_fols

def _get_parser(): 
    parser = ArgumentParser(
        description="""Finds the nearest neighbours for simple structures. 
        Before using on slabs make sure the nn_method works with the bulk 
        structure."""
    )
    parser.add_argument('--hkl', required=True, type=tuple,
    help='Miller index')
    parser.add_argument('-b', '--bulk-energy', required=True, type=str,
    dest='bulk_per_atom', help=('Bulk energy per atom from a converged bulk ' 
    'calculation in eV per atom'))
    parser.add_argument('-p', '--path', default=None, type=str, 
    help='Relative path to the convergence folders (default: cwd)')
    parser.add_argument('--no-enatom', default=True, action='store_false', 
    dest='plt_enatom', help='Turns off energy per atom plotting')
    parser.add_argument('--no-surfen', default=True, action='store_false', 
    dest='plt_surfen', help='Turns off surface energy plotting')
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )

    return parser

def main(): 
    args = _get_parser().parse_args()

    # warnings? 

    parse_fols(args.hkl, args.bulk_per_atom, path=args.path, 
    plt_enatom=args.plt_enatom, plt_surfen=args.plt_surfen, 
    save_csv=args.save_csv)

if __name__ == "__main__":
    main()