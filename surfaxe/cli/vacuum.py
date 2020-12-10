# Misc 
from argparse import ArgumentParser
import os
import warnings 

# Surfaxe 
from surfaxe.vasp_data import vacuum

def _get_parser(): 
    parser = ArgumentParser(
        description="""Gets the energy of the vacuum level."""
    )

    parser.add_argument('-p', '--path', default=None, type=str,  help=('The path ' 
    'to the structure and OUTCAR files (default: cwd)'))

    return parser 

def main(): 
    args = _get_parser().parse_args()
        
    path = os.getcwd()
    if args.path is not None: 
        path = args.path

    vac = vacuum(path)
    print(vac)

if __name__ == "__main__":
    main()