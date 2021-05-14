import unittest
import os
import shutil
from pathlib import Path
from pymatgen.core.surface import Slab
from surfaxe.generation import generate_slabs

ytos = str(Path(__file__).parents[2].joinpath('example_data/generation/CONTCAR_conventional'))
cdte  = str(Path(__file__).parents[2].joinpath('example_data/generation/CdTe.vasp'))

class GenerateSlabsTestCase(unittest.TestCase): 

    def setUp(self): 
        self.ytos = ytos 
        self.cdte = cdte
    
    def test_get_single_hkl(self): 
        ytos_slab = generate_slabs(structure=self.ytos, hkl=(0,0,1), 
        thicknesses=[10,20], vacuums=[10,20], save_slabs=False, max_size=20)

        self.assertEqual(len(ytos_slab), 1)
        self.assertEqual(ytos_slab[0]['slab_index'], 4)
        self.assertWarnsRegex(UserWarning, ('Some generated slabs exceed the '
        'max size specified. Slabs that exceed the max size are: 001_10_10_4'))
        self.assertWarnsRegex(UserWarning, ('Not all combinations of hkl or '
        'slab/vac thicknesses were generated because of repeat structures. '
        'The repeat slabs are: 001_20_10_4, 001_10_20_4, 001_20_20_4'))

    def test_list_hkl(self): 
        ytos_slab = generate_slabs(structure=self.ytos, hkl=[(0,0,1), (1,0,1)], 
        thicknesses=[10,20], vacuums=[10,20], save_slabs=False, max_size=20)

        self.assertEqual(len(ytos_slab), 5)
        self.assertWarnsRegex(UserWarning, 'Not all combinations of hkl or '
        'slab/vac thicknesses were generated because of repeat structures. The '
        'repeat slabs are: 001_10_20_4, 001_20_10_4, 001_20_20_4')
    
    def test_get_max_index(self): 
        ytos_slabs = generate_slabs(structure=self.ytos, 
        hkl=1, thicknesses=[10], vacuums=[10, 20], 
        save_slabs=False, max_size=20)
        
        self.assertIsNotNone(ytos_slabs)
        self.assertEqual(len(ytos_slabs), 7)
        self.assertWarnsRegex(UserWarning, ('Not all combinations of hkl or ' 
        'slab/vac thicknesses were generated because of repeat structures. ' 
        'The repeat slabs are: 001_10_20_4'))
        self.assertWarnsRegex(UserWarning, ('Some generated slabs exceed the ' 
        'max size specified. Slabs that exceed the max size are: 111_10_10_3, '
        '111_10_20_3, 101_10_10_4, 101_10_20_4, 100_10_10_0, 100_10_20_0, '
        '001_10_10_4'))

    def test_get_none(self): 
        ytos_no_slab = generate_slabs(structure=self.ytos, hkl=(0,3,5), 
        thicknesses=[10], vacuums=[10,20], save_slabs=False)

        self.assertEqual(ytos_no_slab, [])
        self.assertWarnsRegex(UserWarning, 'No zero dipole slabs found for '
        'specified Miller index')

    def test_no_structure(self):
        self.assertRaises(FileNotFoundError, generate_slabs, structure='waa', 
        hkl=(0,3,5), thicknesses=[10], vacuums=[10,20], save_slabs=False)
    
    def test_non_centrosymmetric(self): 
        sym_true = generate_slabs(structure=self.cdte, hkl=(1,1,0),
        thicknesses=[10], vacuums=[10,20], save_slabs=False)
        self.assertWarnsRegex(UserWarning, 'Inversion symmetry was not found '
        'in the bulk structure, slabs produced will be non-centrosymmetric')
        
        sym_false = generate_slabs(structure=self.cdte, hkl=(1,1,0),
        thicknesses=[10], vacuums=[10,20], save_slabs=False, is_symmetric=False)

        self.assertEqual(len(sym_true), 2)
        self.assertEqual(len(sym_false), 2)
        self.assertEqual(type(sym_false[0]['slab']), Slab)
    
    def test_save_to_file(self): 
        ytos_slabs = generate_slabs(structure=self.ytos, 
        hkl=(0,0,1), thicknesses=[10], vacuums=[10, 20])
        
        self.assertIsNone(ytos_slabs)
        
        # Check the files created 
        self.assertEqual(len(os.listdir('Y4Ti4S4O10')), 1)
        self.assertIn('POSCAR_001_10_10_15.vasp', os.listdir('Y4Ti4S4O10'))

        # Clean up - get rid of directory created 
        shutil.rmtree('Y4Ti4S4O10')
    
    def test_selective_dynamics(self): 
        ytos_slabs = generate_slabs(structure=self.ytos, 
        hkl=(0,0,1), thicknesses=[30,50], vacuums=[20], layers_to_relax=1)

        self.assertEqual(ytos_slabs[1]['slab'][28].properties['selective_dynamics'], 
        [0.0, 0.0, 0.0])
        self.assertEqual(ytos_slabs[1]['slab'][0].properties['selective_dynamics'], 
        [1.0, 1.0, 1.0])
        self.assertWarnsRegex(UserWarning, 'Some slabs were too thin to fix the '
        'centre of the slab. Slabs with no selective dynamics applied '
        'are: 001_30_20_4')
        self.assertEqual(type(ytos_slabs[0]['slab']), Slab)