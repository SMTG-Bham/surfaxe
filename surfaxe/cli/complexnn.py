# Misc 
from argparse import ArgumentParser
import yaml
import os
import warnings 

# Surfaxe 
from surfaxe.analysis import complex_nn

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
        description="""Finds the nearest neighbours for more complex structures. 
        Uses CutOffDictNN() class as the nearest neighbour method. Check 
        validity on bulk structure before applying to surface slabs."""
    )
    
    parser.add_argument('-s', '--start', default='POSCAR',
    help=('Filename of structure file in any format supported by pymatgen '
          '(default: POSCAR)'))
    parser.add_argument('-a', '--atoms',default=None, nargs='+', type=str,
    help='List of elements in the structure in any order e.g. La Ti O S Ag')
    parser.add_argument('-b', '--bonds', nargs='+',
    dest='cut_off_dict', help='Bond lengths e.g. Bi3+ O2- 2.46 V5+ O2- 1.73')
    parser.add_argument('-e', '--end', default=None,
    help=('Filename of structure file in any format supported by pymatgen. ' 
          'Use if comparing initial and final structures.'))
    parser.add_argument('--oxstates-list', default=None, dest='ox_states_list', 
    help='Add oxidation states to the structure as a list.')
    parser.add_argument('--oxstates-dict', default=None, type=_oxstates_to_dict,
    dest='ox_states_dict', help=('Add oxidation states to the structure as ' 
    'a dictionary e.g. "Fe:3,O:-2"'))
    parser.add_argument('--no-csv', default=True, action='store_false', 
    dest='save_csv', help='Turns off saving data to csv file' )
    parser.add_argument('--csv-fname', default='nn_data.csv',
    dest='csv_fname', help='Filename of the csv file (default: nn_data.csv)')
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

    if args.ox_states_dict: 
        ox_states = args.ox_states_dict 
    elif args.ox_states_list: 
        ox_states = map(float, args.ox_states_list.strip('[]').split(','))
    else: 
        ox_states=None 
    
    # convert cut off dict into correct format 
    cutoff_dict = args.cut_off_dict
    if type(args.cut_off_dict) == list:
        cod = args.cut_off_dict 
        keys, values = ([] for i in range(2))

        for i in [cod[i:i + 3] for i in range(0, len(cod), 3)]: 
            if len(i)==3: 
                values.append(float(i.pop(2)))
            keys.append(tuple(map(str, i)))

        cutoff_dict = dict(zip(keys, values))

    if args.save_csv==True: 
        complex_nn(args.start, args.atoms, cutoff_dict, end=args.end, 
        ox_states=ox_states, save_csv=args.save_csv, csv_fname=args.csv_fname)
    else: 
        nn = complex_nn(args.start, args.atoms, cutoff_dict, end=args.end, 
        ox_states=ox_states, save_csv=args.save_csv, csv_fname=args.csv_fname)
        print(nn)

if __name__ == "__main__":
    main()