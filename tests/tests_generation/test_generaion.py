import unittest
import warnings
import os
from pathlib import Path
from pymatgen import Structure
from pymatgen.core.surface import Slab
from surfaxe.generation import get_all_slabs, get_one_hkl_slabs

ytos = str(Path(__file__).parents[2].joinpath('example_data/generation/CONTCAR_conventional'))
cdte  = str(Path(__file__).parents[2].joinpath('example_data/generation/CdTe.vasp'))

class GetAllTestCase(unittest.TestCase): 

    def setUp(self): 
        self.ytos = ytos
        self.cdte = cdte
    
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

    def test_non_centrosymmetric(self): 
        sym_true = get_all_slabs(structure=self.cdte, max_index=1, 
        thicknesses=[10], vacuums=[10,20], save_slabs=False)
        sym_false = get_all_slabs(structure=self.cdte, max_index=1, 
        thicknesses=[10], vacuums=[10,20], save_slabs=False, is_symmetric=False)

        self.assertEqual(len(sym_true), 0)
        self.assertEqual(len(sym_false), 2)

    def test_save_to_file(self): 
        ytos_slabs = ytos_slabs = get_all_slabs(structure=self.ytos, 
        max_index=1, thicknesses=[10], vacuums=[10, 20])
        
        self.assertIsNone(ytos_slabs)
         
        # add how to actually check if folders are made, check if fmt and poscar 
        # agrs work 

class GetOneTestCase(unittest.TestCase): 

    def setUp(self): 
        self.ytos = ytos 
        self.cdte = cdte
    
    def test_get_one(self): 
        ytos_slab = get_one_hkl_slabs(structure=self.ytos, hkl=(0,0,1), 
        thicknesses=[10,20], vacuums=[10,20], save_slabs=False, max_size=20)

        self.assertEqual(len(ytos_slab), 1)
        self.assertEqual(ytos_slab[0]['s_index'], 4)
        self.assertWarnsRegex(UserWarning, ('Some generated slabs exceed the '
        'max size specified. Slabs that exceed the max size are: 001_10_10_4'))
        self.assertWarnsRegex(UserWarning, ('Not all combinations of hkl or '
        'slab/vac thicknesses were generated because of repeat structures. '
        'The repeat slabs are: 001_20_10_4, 001_10_20_4, 001_20_20_4'))


    def test_get_none(self): 
        ytos_no_slab = get_one_hkl_slabs(structure=self.ytos, hkl=(0,3,5), 
        thicknesses=[10], vacuums=[10,20], save_slabs=False)

        self.assertEqual(ytos_no_slab, [])
    
    #def test_no_structure(self):
       
    #    self.assertRaises(FileNotFoundError, get_one_hkl_slabs(structure='waa', 
    #    hkl=(0,3,5), thicknesses=[10], vacuums=[10,20], save_slabs=False))
    
    def test_non_centrosymmetric(self): 
        sym_true = get_one_hkl_slabs(structure=self.cdte, hkl=(1,1,0),
        thicknesses=[10], vacuums=[10,20], save_slabs=False)
        sym_false = get_one_hkl_slabs(structure=self.cdte, hkl=(1,1,0),
        thicknesses=[10], vacuums=[10,20], save_slabs=False, is_symmetric=False)

        self.assertEqual(len(sym_true), 0)
        self.assertEqual(len(sym_false), 2)
        self.assertEqual(type(sym_false[0]['slab']), Slab)
    
    # add how to actually check if folders are made, check if fmt and poscar 
    # agrs work