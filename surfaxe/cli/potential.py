# Misc 
from argparse import ArgumentParser

from ruamel.yaml.main import YAML

# Surfaxe 
from surfaxe.analysis import electrostatic_potential 

def _get_parser(): 
    parser = ArgumentParser(
        description="""Reads LOCPOT to get the planar and macroscopic 
        potential in specified direction"""
    )

    parser.add_argument('-l', '--locpot', type=str, default='LOCPOT', 
    help='The path to the LOCPOT file (default: ./LOCPOT)')
    parser.add_argument('-p', '--prim-to-conv', type=int, default=1, dest='prim_to_conv',
    help='The number of primitive cells in the conventional cell (default: 1)')
    parser.add_argument('-a', '--axis', type=str, default='c',
    dest='axis', help='Axis of interest; takes abc or xyz (default: c)')
    parser.add_argument('--csv-fname', default='potential.csv', type=str,
    dest='csv_fname', help='Filename of the csv file (default: potential.csv)')
    parser.add_argument('--no-plot', default=True, action='store_false', 
    dest='save_plt', help='Turns off plotting')
    parser.add_argument('-v', '--lattice-vector', type=float, default=None,
    dest='lattice_vector', help='Manually set the periodicity of the slab')
    parser.add_argument('--plt-fname', default='potential.png', type=str,
    dest='plt_fname', help='Filename of the plot (default: potential.png)')
    parser.add_argument('--dpi', default=300, type=int, 
    help='Dots per inch (default: 300)')
    parser.add_argument('-c', '--colors', default=None, nargs='+', type=str, 
    help=('Colors for planar and macroscopic potential plots in any format '
    'supported by mpl e.g. r "#eeefff"; need to supply two valid values where ' 
    'hex colours starting with # need to be surrounded with quotation marks'))
    parser.add_argument('--width', default=6, type=float, 
    help='Width of the figure in inches (default: 6)')
    parser.add_argument('--height', default=5, type=float, 
    help='Height of the figure in inches (default: 5)')
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

        electrostatic_potential(**yaml_args)
        
    else: 
        electrostatic_potential(locpot=args.locpot, prim_to_conv=args.prim_to_conv,
        lattice_vector=args.lattice_vector, axis=args.axis,
        save_csv=True, csv_fname=args.csv_fname, 
        save_plt=args.save_plt, plt_fname=args.plt_fname, dpi=args.dpi, 
        colors=args.colors, width=args.width, height=args.height)


if __name__ == "__main__":
    main()