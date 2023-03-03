# Misc 
from argparse import ArgumentParser
from ruamel.yaml import YAML

# Surfaxe 
from surfaxe.io import plot_bond_analysis

def _get_parser(): 
    parser = ArgumentParser(
        description="""Plots the bond distance with respect to fractional 
        coordinate. Used in conjunction with surfaxe.analysis.bond_analysis."""
    )
    parser.add_argument('-b', '--bond', default=None, nargs='+', type=str,
    help='List of elements e.g. Ti O for a Ti-O bond')
    parser.add_argument('-f', '--filename', default='bond_analysis.csv',
    help='Path to the csv file from bond analysis (default: bond_analysis.csv)')
    parser.add_argument('--plt-fname', default='bond_analysis.png', type=str,
    dest='plt_fname', help='Filename of the plot (default: bond_analysis.png)')
    parser.add_argument('--dpi', default=300, type=int, 
    help='Dots per inch (default: 300)')
    parser.add_argument('-c', '--color', default=None, type=str, 
    help=('Color of the marker in any format supported by mpl e.g. "#eeefff" ' 
    ' hex colours starting with # need to be surrounded with quotation marks' ))
    parser.add_argument('--marker', default='x', type=str, help=('Marker style'))
    parser.add_argument('--markersize', default=5, type=int,help=('Marker size'))
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
        
        plot_bond_analysis(**yaml_args)

    else: 
        plot_bond_analysis(args.bond, filename=args.filename, color=args.color, 
        dpi=args.dpi, width=args.width, height=args.height, 
        plt_fname=args.plt_fname, marker=args.marker, markersize=args.markersize)

if __name__ == "__main__":
    main()