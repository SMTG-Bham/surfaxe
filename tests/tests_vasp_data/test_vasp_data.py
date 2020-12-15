import unittest
import os
from pathlib import Path
from surfaxe.vasp_data import vacuum, core_energy, process_data

data_dir = str(Path(__file__).parents[2].joinpath('example_data/vasp_data'))
analysis_dir = str(Path(__file__).parents[2].joinpath('example_data/analysis'))

class VacuumTestCase(unittest.TestCase): 

    def setUp(self): 
        self.path_csv = os.path.join(data_dir, '101/potential.csv')
        self.path_fol_csv = os.path.join(data_dir, '101')
        self.path_lpt = os.path.join(analysis_dir, 'LOCPOT')
        self.path_fol_lpt = analysis_dir

    def test_vacuum(self): 
        csv = vacuum(self.path_csv)
        fol_csv = vacuum(self.path_fol_csv)
        lpt = vacuum(self.path_lpt)
        fol_lpt = vacuum(self.path_lpt)
        
        self.assertEqual(csv, 7.926)
        self.assertEqual(fol_csv, 7.926)
        self.assertEqual(lpt, 4.557)
        self.assertEqual(fol_lpt, 4.557)
        self.assertWarnsRegex(UserWarning, ('Vacuum electrostatic potential '
        'was not parsed - no LOCPOT or potential.csv files were provided.'), 
        vacuum)
        

class CoreTestCase(unittest.TestCase): 
    def setUp(self): 
        self.outcar = os.path.join(data_dir, '101/OUTCAR')
        self.structure = os.path.join(data_dir, '101/POSCAR')

    def test_core_energy(self): 
        energy = core_energy(core_atom='O', bulk_nn=['Sn', 'Sn', 'Sn'], 
        outcar=self.outcar, structure=self.structure)
        self.assertEqual(energy, -504.2464)

class DataTestCase(unittest.TestCase): 
    def setUp(self): 
        self.path = data_dir 

    def test_data(self): 
        data = process_data(-6.6118, path_to_fols=self.path, save_csv=False) 
        self.assertEqual(data.shape, (2,15))
        self.assertEqual(data['algo'][0], 'Normal')

        data_core = process_data(-6.118, path_to_fols=self.path, save_csv=False, 
        parse_core_energy=True, core_atom='O', bulk_nn=['Sn', 'Sn', 'Sn'])
        self.assertEqual(data_core.shape, (2,16))
        self.assertIn(data_core['core_energy'][0], [-503.5028, -504.2464])