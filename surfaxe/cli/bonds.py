# Misc 
from argparse import ArgumentParser
import yaml
import os
import warnings 

# Surfaxe 
from surfaxe.analysis import bond_analysis
from pymatgen.analysis.local_env import CrystalNN

def _oxstates_to_dict(ox): 
    keys, values = ([] for i in range(2))
    for i in ox.split(','): 
        m = i.split(':')
        values.append(float(m.pop(1)))
        keys.append(str(m.pop(0)))

    ox_states_dict = dict(zip(keys,values))
    return ox_states_dict

def _get_parser(): 
    parser = ArgumentParser(
        description="""Parses the structure looking for bonds between atoms. 
        Check the validity of the nearest neighbour method on the bulk structure 
        before using it on slabs."""
    )

    parser.add_argument('-s', '--structure', default='POSCAR',
    help=('Filename of structure file in any format supported by pymatgen '
          '(default: POSCAR'))
    parser.add_argument('-b', '--bond', default=None, nargs='+', type=str,
    help='List of elements e.g. Ti O for a Ti-O bond')
    parser.add_argument('--oxstates-list', default=None, dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('--oxstates-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', 
    help=('Add oxidation states to the structure as a dictionary '
          ' e.g. "Fe:3,O:-2"'))
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='bond_analysis.csv', 
    dest='csv_fname', help='Filename of the csv file (default: bond_analysis.csv)')
    parser.add_argument('--no-plot', default=True, action='store_false', 
    dest='save_plt', help='Turns off plotting the bond lengths')
    parser.add_argument('--plt-fname', default='bond_analysis.png',
    dest='plt_fname', help='Filename of the plot (default: bond_analysis.png)')
    parser.add_argument('-c', '--color', default=None, type=str, 
    help=('Color of the marker in any format supported by mpl e.g. "#eeefff" ' 
    ' hex colours starting with # need to be surrounded with quotation marks' ))
    parser.add_argument('--width', default=6, type=float, 
    help='Width of the figure in inches (default: 6)')
    parser.add_argument('--height', default=5, type=float, 
    help='Height of the figure in inches (default: 5)')
    parser.add_argument('--dpi', default=300, type=int, 
    help='Dots per inch (default: 300)')
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
 
    if args.ox_states_dict is not None: 
        ox_states = args.ox_states_dict 
    elif args.ox_states_list is not None: 
        ox_states = map(float, args.ox_states_list.strip('[]').split(','))
    else: 
        ox_states=None 

    bond_analysis(args.structure, args.bond, nn_method=CrystalNN(),
    ox_states=ox_states, save_csv=args.save_csv, csv_fname=args.csv_fname, 
    save_plt=args.save_plt, plt_fname=args.plt_fname, dpi=args.dpi, 
    color=args.color, width=args.width, height=args.height)


if __name__ == "__main__":
    main()