import unittest
import warnings
import os
from pathlib import Path
from pymatgen import Structure
from surfaxe.generation import get_all_slabs, get_one_hkl_slabs

ytos = str(Path(__file__).parents[2].joinpath('example_data/generation/CONTCAR_conventional'))

class GetAllTestCase(unittest.TestCase): 

    def setUp(self): 
        self.ytos = ytos
    
    def test_get_all(self): 
        ytos_slabs = get_all_slabs(structure=self.ytos, 
        max_index=1, thicknesses=[10], vacuums=[10, 20], 
        save_slabs=False, max_size=20)
        
        self.assertIsNotNone(ytos_slabs)
        self.assertEqual(len(ytos_slabs), 7)
        self.assertWarnsRegex(UserWarning, ('Not all combinations of hkl or ' 
        'slab/vac thicknesses were generated because of repeat structures. ' 
        'The repeat slabs are: 001_10_20_15'))
        self.assertWarnsRegex(UserWarning, ('Some generated slabs exceed the ' 
        'max size specified. Slabs that exceed the max size are: 111_10_10_3, '
        '101_10_10_9, 100_10_10_10, 001_10_10_15, 111_10_20_3, 101_10_20_9, '
        '100_10_20_10'))

        # add a non-centrosymmetric example test 

class GetOneTestCase(unittest.TestCase): 

    def setUp(self): 
        self.ytos = ytos 
    
    def test_get_one(self): 
        ytos_no_slab = get_one_hkl_slabs(structure=self.ytos, hkl=(0,3,5), 
        thicknesses=[10], vacuums=[10,20], save_slabs=False)

        self.assertEqual(ytos_no_slab, [])

    def test_no_structure(self):
        ytos_no_structure = get_one_hkl_slabs(structure=None, hkl=(0,3,5), 
        thicknesses=[10], vacuums=[10,20], save_slabs=False)
        
        self.assertRaises(FileNotFoundError)