# Misc 
from argparse import ArgumentParser
from ruamel.yaml import YAML

# Surfaxe 
from surfaxe.io import plot_electrostatic_potential

def _get_parser(): 
    parser = ArgumentParser(
        description="""Plots the planar and macroscopic electrostatic 
        potential along one crystallographic direction."""
    )
    parser.add_argument('-f', '--filename', default='potential.csv',
    help='Path to the csv file from bond analysis (default: potential.csv)')
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
        
        plot_electrostatic_potential(**yaml_args)

    else: 
        plot_electrostatic_potential(filename=args.filename, colors=args.colors, 
        dpi=args.dpi, width=args.width, height=args.height, 
        plt_fname=args.plt_fname)

if __name__ == "__main__":
    main()