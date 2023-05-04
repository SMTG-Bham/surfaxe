import os
import unittest
from pathlib import Path
from surfaxe.convergence import parse_energies, parse_structures
import pandas as pd

fols = str(Path(__file__).parents[2].joinpath('example_data/convergence/Y2Ti2S2O5/001'))

class ParseEnergiesTestCase(unittest.TestCase): 

    def setUp(self): 
        self.fols = fols

    def test_parse_energies(self): 
        df = parse_energies(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_surfen=False, save_csv=False)
        b=df[(df['vac_thickness'].astype(int)==30)&(df['slab_thickness'].astype(int)==20)]
        self.assertEqual(df.shape, (6,14))
        self.assertEqual(b['surface_energy_boettger'].values[0], 0.41701338602378163)
        self.assertEqual(b['surface_energy_fm'].values[0], 0.4170133860237656)
        self.assertEqual(b['surface_energy'].values[0], 0.4119406752267468)

    def test_no_pwd(self): 
        self.assertRaises(FileNotFoundError, parse_energies, hkl=(0,0,1), 
        bulk_per_atom=-6.188, save_csv=False, plt_surfen=False)
    
    def test_parse_core_no_atom_set(self): 
        df = parse_energies(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_surfen=False, save_csv=False, 
        parse_core_energy=True)
        self.assertWarnsRegex(UserWarning, 'Core atom or bulk nearest neighbours were not '
            'supplied. Core energy will not be parsed.' )
        self.assertEqual(df.shape, (6,14))

    def test_remove_first_energy(self):
        df = parse_energies(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_surfen=False, save_csv=False, 
        remove_first_energy=True)
        self.assertWarnsRegex(UserWarning, 'First data point was not removed '
        '- less than three data points were present in dataset')
    
    def test_parse_energie_no_parallel(self): 
        df = parse_energies(hkl=(0,0,1), bulk_per_atom=-8.83099767, 
        path_to_fols=self.fols, plt_surfen=False, save_csv=False, processes=1)
        b=df[(df['vac_thickness'].astype(int)==30)&(df['slab_thickness'].astype(int)==20)]
        self.assertEqual(df.shape, (6,14))
        self.assertEqual(b['surface_energy_boettger'].values[0], 0.41701338602378163)
        self.assertEqual(b['surface_energy_fm'].values[0], 0.4170133860237656)
        self.assertEqual(b['surface_energy'].values[0], 0.4119406752267468)

class ParseStructuresTestCase(unittest.TestCase):

    def setUp(self): 
        self.fols = fols

    def tearDown(self) -> None:
        if os.path.isfile('Y2Ti2S2O5_parsed_metadata.json'): 
            os.remove('Y2Ti2S2O5_parsed_metadata.json')
        return super().tearDown()

    def test_parse_structures_json(self): 
        parse_structures((0,0,1), structure_file='POSCAR', path_to_fols=self.fols)

        self.assertIn('Y2Ti2S2O5_parsed_metadata.json', os.listdir(os.getcwd()))

    def test_parse_bond(self): 
        d = parse_structures((0,0,1), structure_file='POSCAR', path_to_fols=self.fols, save_json=False,
        bond=['Y', 'O'])

        self.assertIn('Y-O', d[0]['bonds'].keys())
        self.assertEqual(d[0]['bonds']['Y-O']['Y_c_coord'][0], 0.708615)
        self.assertEqual(d[0]['bonds']['Y-O']['Y-O_bond_distance'][3], 2.4095195472222635)    

    def test_parse_bonds(self): 
        d = parse_structures((0,0,1), structure_file='POSCAR', path_to_fols=self.fols, save_json=False,
        bond=[['Y', 'O'], ['Ti', 'O']])

        self.assertIn('Y-O', d[0]['bonds'].keys())
        self.assertIn('Ti-O', d[0]['bonds'].keys())
        
        self.assertEqual(d[0]['bonds']['Y-O']['Y_c_coord'][0], 0.708615)
        self.assertEqual(d[0]['bonds']['Y-O']['Y-O_bond_distance'][3], 2.4095195472222635)

        self.assertEqual(len(d[0]['bonds']['Ti-O']['Ti_c_coord']), 4)
        self.assertAlmostEqual(d[0]['bonds']['Ti-O']['Ti-O_bond_distance'][1], 1.9193156521064114)
        self.assertEqual(d[0]['bonds']['Ti-O']['Ti_c_coord'][3], 0.415056)

    def test_parse_bonds_auto(self): 
        d = parse_structures((0,0,1), structure_file='POSCAR', path_to_fols=self.fols, save_json=False,
        bond='auto')

        self.assertIn('Y-O', d[0]['bonds'].keys())
        self.assertIn('Ti-O', d[0]['bonds'].keys())
        self.assertIn('Ti-S', d[0]['bonds'].keys())
        self.assertIn('Y-S', d[0]['bonds'].keys())
        
        self.assertEqual(d[0]['bonds']['Y-O']['Y_c_coord'][0], 0.708615)
        self.assertEqual(d[0]['bonds']['Y-O']['Y-O_bond_distance'][3], 2.4095195472222635)
  
        self.assertEqual(len(d[0]['bonds']['Ti-O']['Ti_c_coord']), 4)
        self.assertAlmostEqual(d[0]['bonds']['Ti-O']['Ti-O_bond_distance'][1], 1.9193156521064114)
        self.assertEqual(d[0]['bonds']['Ti-O']['Ti_c_coord'][3], 0.415056)


    def test_parse_bonds_none(self): 
        d = parse_structures((0,0,1), structure_file='POSCAR', path_to_fols=self.fols, save_json=False,
        bond=None)

        self.assertNotIn('bonds', d[0].keys())
    
    def test_parse_bonds_nonsense(self): 
        with self.assertRaises(TypeError) as e: 
            parse_structures((0,0,1), structure_file='POSCAR', path_to_fols=self.fols, save_json=False,
            bond='nonsense')

    
    # should also include separate tests for parse_vacuum but the problem is the 
    # locpots needed for that are large even for small systems 