# Misc 
from argparse import ArgumentParser
import json
import os
import warnings 

# Surfaxe 
from surfaxe.generation import get_all_slabs


def _get_parser(): 
    parser = ArgumentParser(
        description="""Generates all unique slabs with specified maximum 
        Miller index, minimum slab and vacuum thicknesses. It includes all 
        combinations for multiple zero dipole symmetric terminations for the 
        same Miller index"""
    )

    parser.add_argument('-s', '--structure', required=True, type=str,
    help='Filename of structure file in any format supported by pymatgen')
    parser.add_argument('--hkl', required=True, type=int,
    help='The maximum Miller index to go up to.')
    parser.add_argument('-t', '-thicknesses', required=True, type=list,
    help='The minimum size of the slab in Angstroms.')
    parser.add_argument('-v', '--vacuums', required=True, type=list,
    help='The minimum size of the vacuum in Angstroms.')
    parser.add_argument('-r', '--fols', default=False, action='store_true', 
    help=('Makes folders for each termination and slab/vacuum thickness ' 
          'combinations containing POSCARs (default: False)'))
    parser.add_argument('-f', '--files', default=False, action='store_true', 
    help='Makes INCAR, POTCAR and KPOINTS files in each folder (default: False)')
    parser.add_argument('--max-size', default=500, dest='max_size', type=int,
    help=('The maximum number of atoms in the slab specified to raise warning ' 
          'about slab size. Even if the warning is raised, it still outputs ' 
          'the slabs regardless.'))
    parser.add_argument('-b', '--bonds', default=None, type=dict,
    help='Bonds to keep intact while cleaving the slab')
    parser.add_argument('-c', '--center-slab', default=True, dest='center_slab',
    action='store_false', help='The position of the slab in the simulation cell')
    parser.add_argument('--oxstates', default=None, dest='ox_states',
    help='Add oxidation states to the structure.')
    parser.add_argument('--no-save', default=True, action='store_false', 
    dest='save_slabs', 
    help='Whether to save the slabs to file (default: True)')
    parser.add_argument('--no-sym', default=True, action='store_false', 
    dest='sym', help='Whether the slabs cleaved should have inversion symmetry.')
    parser.add_argument('--config-dict', default='PBEsol_config.json', dest='config_dict',
    help='Specifies the dictionary used for the generation of the input files')
    parser.add_argument('-i', '--incar', default=None,
    help='Overrides the default INCAR parameter settings')
    parser.add_argument('-k', '--kpoints', default=None,
    help='Overrides the default KPOINTS settings.')
    parser.add_argument('-p', '--potcar', default=None,
    help='Overrides the default POTCAR settings')

    return parser

def main(): 
    args = _get_parser().parse_args()

    # Warnings? 

    get_all_slabs(args.structure, args.hkl, args.thicknesses, args.vacuums, 
    make_fols=args.fols, make_input_files=args.files, max_size=args.max_size, 
    bonds=args.bonds, center_slab=args.center_slab, ox_states=args.ox_states, 
    save_slabs=args.save_slabs, is_symmetric=args.sym, config_dict=args.config_dict, 
    user_incar_settings=args.incar, user_potcar_settings=args.potcar, 
    user_kpoints_settings=args.kpoints)

if __name__ == "__main__":
    main()