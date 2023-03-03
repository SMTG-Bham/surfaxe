# Misc 
from argparse import ArgumentParser
from ruamel.yaml import YAML

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
    parser.add_argument('--oxi-list', default=None, dest='ox_states_list', 
    nargs='+', type=float, 
    help='Add oxidation states to the structure as a list e.g. 3 3 -2 -2 -2')
    parser.add_argument('--oxi-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', help=('Add oxidation states to the structure as ' 
    'a dictionary e.g. Fe:3,O:-2'))
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Prints data to terminal' )
    parser.add_argument('--csv-fname', default='bond_analysis.csv', 
    dest='csv_fname', help='Filename of the csv file (default: bond_analysis.csv)')
    parser.add_argument('--no-plot', default=True, action='store_false', 
    dest='save_plt', help='Turns off plotting')
    parser.add_argument('--plt-fname', default='bond_analysis.png',
    dest='plt_fname', help='Filename of the plot (default: bond_analysis.png)')
    parser.add_argument('-c', '--color', default=None, type=str, 
    help=('Color of the marker in any format supported by mpl e.g. "#eeefff" ' 
    ' hex colours starting with # need to be surrounded with quotation marks' ))
    parser.add_argument('--marker', default='x', type=str, help=('Marker style'))
    parser.add_argument('--markersize', default=5, type=int,help=('Marker size'))
    parser.add_argument('--width', default=6, type=float, 
    help='Width of the figure in inches (default: 6)')
    parser.add_argument('--height', default=5, type=float, 
    help='Height of the figure in inches (default: 5)')
    parser.add_argument('--dpi', default=300, type=int, 
    help='Dots per inch (default: 300)')
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
                
        ba = bond_analysis(**yaml_args)
        if ('save_csv', False) in yaml_args.items(): 
            print(ba)
    
    else: 
        if args.ox_states_dict is not None: 
            ox_states = args.ox_states_dict 
        elif args.ox_states_list is not None: 
            ox_states = args.ox_states_list
        else: 
            ox_states=None 
        

        ba = bond_analysis(args.structure, args.bond, nn_method=CrystalNN(),
        ox_states=ox_states, save_csv=args.save_csv, csv_fname=args.csv_fname, 
        save_plt=args.save_plt, plt_fname=args.plt_fname, dpi=args.dpi, 
        color=args.color, width=args.width, height=args.height, 
        marker=args.marker, markersize=args.markersize)
        
        if not args.save_csv: 
            print(ba)


if __name__ == "__main__":
    main()