# Misc 
from argparse import ArgumentParser
import os

from ruamel.yaml.main import YAML

# Surfaxe 
from surfaxe.convergence import parse_structures

def _get_parser(): 
    parser = ArgumentParser(
        description="""Parses the convergence folders to save the structures to 
        json and optionally does bond analysis for one bond. Bond analysis 
        works only for systems where CrystalNN can reliably determine the 
        coordination environment. 
       """
    )
    parser.add_argument('--hkl', help='Miller index e.g. 0,0,-1')
    parser.add_argument('-s', '--structure', default=None, type=str,
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('-b', '--bond', default=None, nargs='+', type=str,
    help='List of elements e.g. Ti O for a Ti-O bond')
    parser.add_argument('-p', '--path', default=None, type=str, 
    help='Relative path to the convergence folders (default: cwd)')
    parser.add_argument('--json-fname', default=None, type=str,
    dest='json_fname', help=('Filename of the json file (default: '
    'formula_parsed_metadata.json)'))
    parser.add_argument('-v', '--verbose', default=False, action='store_true', 
    help=('Whether or not to print extra info about the folders being parsed.'
    ' (default: False)'))
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
        
        parse_structures(**yaml_args)

    else: 
        if not args.hkl: 
            raise ValueError('hkl was not supplied')

        hkl = tuple(map(int, args.hkl.strip('[]()').split(',')))
        
        path = os.getcwd()
        if args.path is not None: 
            path = args.path

        parse_structures(hkl, args.bulk_per_atom, path_to_fols=path, 
        parse_core_energy=args.parse_core, core_atom=args.core, bulk_nn=args.nn, 
        parse_vacuum=args.parse_vacuum, plt_enatom=args.plt_enatom, 
        plt_surfen=args.plt_surfen, csv_fname=args.json_fname, 
        verbose=args.verbose)

if __name__ == "__main__":
    main()