import unittest
import warnings
import os
from pymatgen import Structure
from surfaxe.data import vacuum, core_energy

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(THIS_DIR, os.pardir, 'data/')

class VacuumTestCase(unittest.TestCase): 

    def setUp(self): 
        self.path = data_dir

    def test_vacuum(self): 
        value = vacuum(self.path)
        self.assertEqual(value, 4.557)        

class CoreTestCase(unittest.TestCase): 
    def setUp(self): 
        self.path = data_dir

    def test_core_energy(self): 
        energy = core_energy(self.path, 'O', ['Ti', 'Ti', 'Y', 'Y'])
        self.assertEqual(energy, -505.7559)

class DataTestCase(unittest.TestCase): 
    def setUp(self): 
        pass 
    def test_data(self): 
        pass 