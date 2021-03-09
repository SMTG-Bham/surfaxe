# Misc 
from argparse import ArgumentParser
import yaml
import os
import warnings 

# Surfaxe 
from surfaxe.analysis import electrostatic_potential 

def _get_parser(): 
    parser = ArgumentParser(
        description="""Reads LOCPOT to get the planar and macroscopic 
        potential in specified direction"""
    )

    parser.add_argument('-v', '--lattice-vector', type=float, 
    dest='lattice_vector', help='The periodicity of the slab')
    parser.add_argument('-l', '--locpot', type=str, default='LOCPOT', 
    help='The path to the LOCPOT file (default: ./LOCPOT)')
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='potential.csv', type=str,
    dest='csv_fname', help='Filename of the csv file (default: potential.csv)')
    parser.add_argument('--no-plot', default=True, action='store_false', 
    dest='save_plt', help='Turns off plotting')
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
            yaml_args = yaml.safe_load(y)

        ep = electrostatic_potential(**yaml_args)
        if ('save_csv', False) in yaml_args.items(): 
            print(ep)
        
    else: 
        ep = electrostatic_potential(args.lattice_vector, locpot=args.locpot, 
        axis=2, save_csv=args.save_csv, csv_fname=args.csv_fname, 
        save_plt=args.save_plt, plt_fname=args.plt_fname, dpi=args.dpi, 
        colors=args.colors, width=args.width, height=args.height)

        if not args.save_csv: 
            print(ep)

if __name__ == "__main__":
    main()