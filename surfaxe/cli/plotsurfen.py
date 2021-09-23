# Misc 
from argparse import ArgumentParser
import pandas as pd
from ruamel.yaml.main import YAML

# Surfaxe 
from surfaxe.io import plot_surfen

def _get_parser(): 
    parser = ArgumentParser(
        description="""Plots the surface energy for all terminations."""
    )
    parser.add_argument('-f', '--filename', 
    help='Path to the csv file from parsefols with data')
    parser.add_argument('--plt-fname', default='surface_energy.png', type=str,
    dest='plt_fname', help='Filename of the plot (default: surface_energy.png)')
    parser.add_argument('--dpi', default=300, type=int, 
    help='Dots per inch (default: 300)')
    parser.add_argument('-c', '--colors', default=None, nargs='+', type=str, 
    help=('Colours for different vacuum thicknesses plots in any format '
    'supported by mpl e.g. r g "#eeefff" where hex colours starting with # need '
    'to be surrounded with quotation marks' ))
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
        
        df = pd.read_csv(yaml_args['filename'])
        plot_surfen(df=df, **yaml_args)
    
    else: 
        df = pd.read_csv(args.filename)
        plot_surfen(df, colors=args.colors, dpi=args.dpi, width=args.width, 
        height=args.height,  plt_fname=args.plt_fname)

if __name__ == "__main__":
    main()