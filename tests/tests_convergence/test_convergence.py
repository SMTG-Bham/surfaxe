import unittest
import warnings
import os
from pathlib import Path
from pymatgen import Structure
from pymatgen.core.surface import Slab
from surfaxe.convergence import parse_fols

fols = str(Path(__file__).parents[2].joinpath('example_data/convergence/Y2Ti2S2O5'))

class ParseFolsTestCase(unittest.TestCase): 

    def setUp(self): 
        self.fols = fols

    def test_parse_fols(self): 
        parse_fols_data = parse_fols(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_enatom=False, plt_surfen=False, 
        save_csv=False)

        self.assertEqual(parse_fols_data.shape, (6,12))

    def test_no_pwd(self): 
        self.assertRaises(FileNotFoundError, parse_fols, hkl=(0,0,1), 
        bulk_per_atom=-6.188, save_csv=False, plt_enatom=False, plt_surfen=False)
        