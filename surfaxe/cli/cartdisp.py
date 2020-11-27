# Misc 
from argparse import ArgumentParser
import json
import os
import warnings 

# Surfaxe 
from surfaxe.analysis import cart_displacements

def _get_parser(): 
    parser = ArgumentParser(
        description="""Finds the nearest neighbours for simple structures. 
        Before using on slabs make sure the nn_method works with the bulk 
        structure."""
    )

    parser.add_argument('-s' '--start', required=True, type=str, 
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('-e', '--end', required=True, type=str,
    help=('Filename of structure file in any format supported by pymatgen. ' 
          'Use if comparing initial and final structures.'))
    parser.add_argument('-a', '--atoms', required=True, type=list, dest='elements',
    help='List of elements in the structure in any order')
    parser.add_argument('--max-disp', type=float, default=0.1, dest='max_disp', 
    help='The maximum displacement shown (default: 0.1')
    parser.add_argument('--no-txt', default=True, action='store_false', 
    dest='save_txt', help='Turns off saving data to a txt file' )
    parser.add_argument('--txt-fname', default='cart_displacament.txt', type=str,
    dest='txt_fname', help=('Filename of the txt file (default: '
    'cart_displacement.txt)'))

    return parser 

def main(): 
    args = _get_parser().parse_args()

    # warnings? 
    cart_displacements(args.start, args.end, args.elements, 
    max_disp=args.max_disp, save_txt=args.save_txt, txt_fname=args.txt_fname)

if __name__ == "__main__":
    main()