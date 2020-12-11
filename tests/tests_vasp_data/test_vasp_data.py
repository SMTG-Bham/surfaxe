import unittest
import warnings
import os
from pathlib import Path
from pymatgen import Structure
from surfaxe.vasp_data import vacuum, core_energy, process_data

data_dir = str(Path(__file__).parents[1].joinpath('example_data/vasp_data'))

class VacuumTestCase(unittest.TestCase): 

    def setUp(self): 
        self.path = os.path.join(data_dir, '101/potential.csv')

    def test_vacuum(self): 
        value = vacuum(self.path)
        self.assertEqual(value, 7.926)        

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
        pass 
    def test_data(self): 
        pass 