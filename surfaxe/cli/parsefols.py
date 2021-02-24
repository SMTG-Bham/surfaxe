# Misc 
from argparse import ArgumentParser
import yaml
import os
import warnings 

# Surfaxe 
from surfaxe.convergence import parse_fols

def _get_parser(): 
    parser = ArgumentParser(
        description="""Parses the convergence folders to get the surface energy, 
        total energy, energy per atom, band gap and time taken for each slab and 
        vacuum thickness combination. It can optionally parse vacuum and core 
        level energies. Core energy parsing works only for systems where 
        CrystalNN can reliably determine the coordination environment. 
        By default convergence plots are turned off - they can be customised by 
        using surfaxe-plot-enatom and surfaxe-plot-surfen"""
    )
    parser.add_argument('--hkl', help='Miller index e.g. 0,0,-1')
    parser.add_argument('-b', '--bulk-energy', type=float,
    dest='bulk_per_atom', help=('Bulk energy per atom from a converged bulk ' 
    'calculation in eV per atom'))
    parser.add_argument('-p', '--path', default=None, type=str, 
    help='Relative path to the convergence folders (default: cwd)')
    parser.add_argument('--parse-core', default=False, action='store_true', 
    dest='parse_core', help=('Attempts to parse core energies from a supplied '
    'OUTCAR (default: False)'))
    parser.add_argument('--core', default=None, type=str, 
    help='The symbol of atom the core state energy should be parsed from')
    parser.add_argument('--nn', default=None, nargs='+', type=str,
    help='The symbols of the core atom nearest neighbours e.g. Y Y Ti Ti' )
    parser.add_argument('--parse-vacuum', default=False, action='store_true', 
    dest='parse_vacuum', help=('Attempts to get the maximum value of planar '
    'potential from a LOCPOT (default: False)'))
    parser.add_argument('--plot-enatom', default=False, action='store_true', 
    dest='plt_enatom', 
    help='Plot basic energy per atom vs slab thickness figure (default: False)')
    parser.add_argument('--plot-surfen', default=False, action='store_true', 
    dest='plt_surfen', 
    help='Plot basic surface energy vs slab thickness figure (default: False)')
    parser.add_argument('--csv-fname', default=None, type=str,
    dest='csv_fname', help='Filename of the csv file (default: hkl_data.csv)')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', 
    help=('Whether or not to print extra info about the folders being parsed.'
    ' (default: False)'))
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

    if not args.hkl or not args.bulk_per_atom: 
        raise ValueError('hkl or bulk energy per atom were not supplied')

    path = os.getcwd()
    if args.path is not None: 
        path = args.path

    parse_fols(hkl, args.bulk_per_atom, path_to_fols=path, 
    parse_core_energy=args.parse_core, core_atom=args.core, bulk_nn=args.nn, 
    parse_vacuum=args.parse_vacuum, plt_enatom=args.plt_enatom, 
    plt_surfen=args.plt_surfen, save_csv=True, csv_fname=args.csv_fname, 
    verbose=args.verbose)

if __name__ == "__main__":
    main()