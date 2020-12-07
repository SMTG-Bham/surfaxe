# Misc 
from argparse import ArgumentParser
import yaml
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
    parser.add_argument('--hkl', help='Miller index e.g. 0,0,-1')
    parser.add_argument('-b', '--bulk-energy', type=float,
    dest='bulk_per_atom', help=('Bulk energy per atom from a converged bulk ' 
    'calculation in eV per atom'))
    parser.add_argument('-p', '--path', default=None, type=str, 
    help='Relative path to the convergence folders (default: cwd)')
    parser.add_argument('--no-enatom', default=True, action='store_false', 
    dest='plt_enatom', help='Turns off energy per atom plotting')
    parser.add_argument('--no-surfen', default=True, action='store_false', 
    dest='plt_surfen', help='Turns off surface energy plotting')
    parser.add_argument('--yaml', default=False, action='store_true', 
    help='Read optional args from surfaxe_config.yaml file.')

    return parser

def main(): 
    args = _get_parser().parse_args()
    if args.hkl is not None: 
        hkl = tuple(map(int, args.hkl.strip('[]()').split(',')))

    if args.yaml==True: 
        with open('surfaxe_config.yaml', 'r') as y: 
            yaml_args = yaml.load(y)
        args.update(
            (k, yaml_args[k]) for k in args.keys() and yaml_args.keys()
        )

    path = os.getcwd()
    if args.path is not None: 
        path = args.path

    parse_fols(hkl, args.bulk_per_atom, path=path, 
    plt_enatom=args.plt_enatom, plt_surfen=args.plt_surfen, 
    save_csv=True)

if __name__ == "__main__":
    main()