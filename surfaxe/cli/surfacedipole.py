from ruamel.yaml.main import YAML
from argparse import ArgumentParser
# Surfaxe 
from surfaxe.analysis import surface_dipole 

def _get_parser(): 
    parser = ArgumentParser(
        description="""Calculates the surface dipole moment of a slab for band alignments"""
    )

    parser.add_argument('-f', '--filename', type=str, default='LOCPOT', 
    help='The path to the LOCPOT or parsed csv file (default: ./LOCPOT)')
    parser.add_argument('-p', '--prim-to-conv', type=int, default=1, dest='prim_to_conv',
    help='The number of primitive cells in the conventional cell (default: 1)')
    parser.add_argument('-a', '--axis', type=str, default='c',
    dest='axis', help='Axis of interest; takes abc or xyz (default: c)')
    parser.add_argument('-v', '--lattice-vector', type=float, default=None,
    dest='lattice_vector', help='Manually set the periodicity of the slab')
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

        ep = surface_dipole(**yaml_args)
        print(ep)
        
    else: 
        ep = surface_dipole(filename=args.filename, prim_to_conv=args.prim_to_conv, lattice_vector=args.lattice_vector, axis=args.axis,
        save_csv=False, save_plt=False, )
        print(ep)

if __name__ == "__main__":
    main()