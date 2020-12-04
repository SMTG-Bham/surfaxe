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
    parser.add_argument('--hkl', required=True, help='Miller index')
    parser.add_argument('-b', '--bulk-energy', required=True, type=float,
    dest='bulk_per_atom', help=('Bulk energy per atom from a converged bulk ' 
    'calculation in eV per atom'))
    parser.add_argument('-p', '--path', default=None, type=str, 
    help='Relative path to the convergence folders (default: cwd)')
    parser.add_argument('--no-enatom', default=True, action='store_false', 
    dest='plt_enatom', help='Turns off energy per atom plotting')
    parser.add_argument('--no-surfen', default=True, action='store_false', 
    dest='plt_surfen', help='Turns off surface energy plotting')

    return parser

def main(): 
    args = _get_parser().parse_args()
    hkl = tuple(map(int, args.hkl.strip('[]()').split(',')))

    path = os.getcwd()
    if args.path is not None: 
        path = args.path

    parse_fols(hkl, args.bulk_per_atom, path=path, 
    plt_enatom=args.plt_enatom, plt_surfen=args.plt_surfen, 
    save_csv=True)

if __name__ == "__main__":
    main()