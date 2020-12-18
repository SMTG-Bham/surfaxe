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
    help='The path to the LOCPOT file (default: ./LOCPOT ')
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='potential.csv', type=str,
    dest='csv_fname', help='Filename of the csv file')
    parser.add_argument('--no-plot', default=True, action='store_false', 
    dest='save_plt', help='Turns off plotting the potentials')
    parser.add_argument('--plt-fname', default='potential.png', type=str,
    dest='plt_fname', help='Filename of the plot')
    parser.add_argument('--dpi', default=300, type=int, help='Dots per inch')
    parser.add_argument('--yaml', default=False, action='store_true', 
    help='Read optional args from surfaxe_config.yaml file.')

    return parser

def main(): 
    args = _get_parser().parse_args()

    if args.yaml==True: 
        with open('surfaxe_config.yaml', 'r') as y: 
            yaml_args = yaml.load(y)
        args.update(
            (k, yaml_args[k]) for k in args.keys() and yaml_args.keys()
        ) 

    electrostatic_potential(args.lattice_vector, locpot=args.locpot, 
    axis=2, save_csv=args.save_csv, csv_fname=args.csv_fname, 
    save_plt=args.save_plt, plt_fname=args.plt_fname, dpi=args.dpi)

if __name__ == "__main__":
    main()