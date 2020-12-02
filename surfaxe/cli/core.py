# Misc 
from argparse import ArgumentParser
import os
import warnings 

# Surfaxe 
from surfaxe.data import core

def _get_parser(): 
    parser = ArgumentParser(
        description="""Parses the structure and OUTCAR files for the core level 
        energy. Check the validity of nearest neighbour method on the bulk 
        structure before using it on slabs."""
    )

    parser.add_argument('-p', '--path', default=None, type=str,  help=('The path ' 
    'to the structure and OUTCAR files'))
    parser.add_argument('-a', '--atom', required=True, type=str, dest='core_atom', 
    help='The symbol of atom the core state energy level should be parsed from')
    parser.add_argument('-b', '--bulknn', required=True, type=list, 
    dest='bulk_nn', help='The symbols of the nearest neighbours of the core atom')
    parser.add_argument('-o', '--orbital', type=str, default='1s', 
    help='The orbital of core state (default: 1s)')
    # not sure how to pass class as the default here? 
    parser.add_argument('-n', '--nnmethod', default='CrystalNN()', type=str, 
    dest='nn_method', 
    help='The pymatgen local_env nearest neighbour method (default: CrystalNN()')
    parser.add_argument('--oxstates', default=None, dest='ox_states',
    help='Add oxidation states to the structure.')
    parser.add_argument('-s', '--structure', default='vasprun.xml', type=str,
    help=('Filename of structure file in any format supported by pymatgen '
    'default: vasprun.xml'))

    return parser

def main(): 
    args = _get_parser().parse_args()

    path = os.getcwd()
    if args.path is not None: 
        path = args.path

    core(path, args.core_atom, args.bulk_nn, 
    orbital=args.orbital, ox_states=args.ox_states, nn_method=args.nn_method, 
    structure=args.structure)

if __name__ == "__main__":
    main()