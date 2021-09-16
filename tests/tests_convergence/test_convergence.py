import unittest
from pathlib import Path
from surfaxe.convergence import parse_energies
import pandas as pd

fols = str(Path(__file__).parents[2].joinpath('example_data/convergence/Y2Ti2S2O5/001'))

class ParseEnergiesTestCase(unittest.TestCase): 

    def setUp(self): 
        self.fols = fols

    def test_parse_energies(self): 
        df = parse_energies(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_surfen=False, save_csv=False)
        b=df[(df['vac_thickness'].astype(int)==30)&(df['slab_thickness'].astype(int)==20)]
        self.assertEqual(df.shape, (6,15))
        self.assertEqual(b['surface_energy_boettger'][0], 0.41701338602378163)
        self.assertEqual(b['surface_energy_fm'][0], 0.4170133860237656)
        self.assertEqual(b['surface_energy'][0], 0.4119406752267468)

    def test_no_pwd(self): 
        self.assertRaises(FileNotFoundError, parse_energies, hkl=(0,0,1), 
        bulk_per_atom=-6.188, save_csv=False, plt_surfen=False)
    
    def test_parse_core_no_atom_set(self): 
        df = parse_energies(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_surfen=False, save_csv=False, 
        parse_core_energy=True)
        self.assertWarnsRegex(UserWarning, 'Core atom or bulk nearest neighbours were not '
            'supplied. Core energy will not be parsed.' )
        self.assertEqual(df.shape, (6,15))

    def test_remove_first_energy(self):
        df = parse_energies(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_surfen=False, save_csv=False, 
        remove_first_energy=True)
        self.assertWarnsRegex(UserWarning, 'First data point was not removed '
        '- less than three data points were present in dataset')