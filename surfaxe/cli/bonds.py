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

def _bonds_to_list(bonds): 
    bonds_list = []
    for bond in bonds.split(','): 
        waa = []
        for el in bond.split('.'): 
            waa.append(el)
        bonds_list.append(waa)

    return bonds_list 

def _get_parser(): 
    parser = ArgumentParser(
        description="""Parses the structure looking for bonds between atoms. 
        Check the validity of the nearest neighbour method on the bulk structure 
        before using it on slabs."""
    )

    parser.add_argument('-s', '--structure',
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('-b', '--bonds', type=_bonds_to_list, dest='list_of_bonds', 
    help='List of bonds as lists to compare in any order (e.g. Y-O,Ti-S')
    parser.add_argument('--oxstates-list', default=None, dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('--oxstates-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', help=('Add oxidation states to the structure as ' 
    'a dictionary e.g. "Fe:3,O:-2"'))
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='bond_analysis.csv', type=str,
    dest='csv_fname', help='Filename of the csv file (default: bond_analysis.csv)')
    parser.add_argument('--no-plot', default=True, action='store_false', 
    dest='save_plt', help='Turns off plotting the bond lengths')
    parser.add_argument('--plt-fname', default='bond_analysis.png', type=str,
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
 
    if args.ox_states_dict is not None: 
        ox_states = args.ox_states_dict 
    elif args.ox_states_list is not None: 
        ox_states = map(float, args.ox_states_list.strip('[]').split(','))
    else: 
        ox_states=None 

    bond_analysis(args.structure, args.list_of_bonds, nn_method=CrystalNN(),
    ox_states=ox_states, save_csv=args.save_csv, csv_fname=args.csv_fname, 
    save_plt=args.save_plt, plt_fname=args.plt_fname, dpi=args.dpi)


if __name__ == "__main__":
    main()