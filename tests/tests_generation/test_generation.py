import unittest
import os
import shutil
from pathlib import Path
from pymatgen.core.surface import Slab
from pymatgen.core import Structure
from surfaxe.generation import generate_slabs

ytos = str(Path(__file__).parents[2].joinpath('example_data/generation/CONTCAR_conventional'))
cdte  = str(Path(__file__).parents[2].joinpath('example_data/generation/CdTe.vasp'))

class GenerateSlabsTestCase(unittest.TestCase): 

    def setUp(self): 
        self.ytos = ytos 
        self.cdte = cdte
        self.ytos_pmg = Structure.from_file(ytos)
    
    def test_get_single_hkl(self): 
        ytos_slab = generate_slabs(structure=self.ytos, hkl=(0,0,1), 
        thicknesses=[10,20], vacuums=[10,20], save_slabs=False, 
        save_metadata=False, max_size=20)

        self.assertEqual(len(ytos_slab), 1)
        self.assertEqual(ytos_slab[0]['slab_index'], 4)
        self.assertWarnsRegex(UserWarning, ('Some generated slabs exceed the '
        'max size specified. Slabs that exceed the max size are: 001_10_10_4'))
        self.assertWarnsRegex(UserWarning, ('Not all combinations of hkl or '
        'slab/vac thicknesses were generated because of repeat structures. '
        'The repeat slabs are: 001_20_10_4, 001_10_20_4, 001_20_20_4'))

    def test_list_hkl(self): 
        ytos_slab = generate_slabs(structure=self.ytos, hkl=[(0,0,1), (1,0,1)], 
        thicknesses=[10,20], vacuums=[10,20], save_slabs=False, max_size=20, 
        save_metadata=False)

        self.assertEqual(len(ytos_slab), 5)
        self.assertWarnsRegex(UserWarning, 'Not all combinations of hkl or '
        'slab/vac thicknesses were generated because of repeat structures. The '
        'repeat slabs are: 001_10_20_4, 001_20_10_4, 001_20_20_4')
    
    def test_get_max_index(self): 
        ytos_slabs = generate_slabs(structure=self.ytos, 
        hkl=1, thicknesses=[10], vacuums=[10, 20], 
        save_slabs=False, max_size=20, save_metadata=False)
        
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
        with self.assertRaises(ValueError) as cm: 
            slabs = generate_slabs(structure=self.ytos, hkl=(0,3,5), 
        thicknesses=[10], vacuums=[10,20], save_slabs=False, save_metadata=False)
        
        ex = cm.exception
        self.assertEqual(str(ex),'No zero dipole (Tasker I or II) slabs found for specified '
        'Miller index')
            


    def test_no_structure(self):
        self.assertRaises(FileNotFoundError, generate_slabs, structure='waa', 
        hkl=(0,3,5), thicknesses=[10], vacuums=[10,20], save_slabs=False)

        self.assertRaises(TypeError, generate_slabs, structure=1, 
        hkl=(0,3,5), thicknesses=[10], vacuums=[10,20], save_slabs=False)


    def test_pmg_structure(self): 
        pmg_struc = generate_slabs(structure=self.ytos_pmg, 
        hkl=1, thicknesses=[10], vacuums=[10, 20], 
        save_slabs=False, save_metadata=False)

        self.assertIsNotNone(pmg_struc)
        self.assertEqual(len(pmg_struc), 7)
    
    def test_non_centrosymmetric(self): 
        sym_true = generate_slabs(structure=self.cdte, hkl=(1,1,0),
        thicknesses=[10], vacuums=[10,20], save_slabs=False, save_metadata=False)
        self.assertWarnsRegex(UserWarning, 'Inversion symmetry was not found '
        'in the bulk structure, slabs produced will be non-centrosymmetric')
        
        sym_false = generate_slabs(structure=self.cdte, hkl=(1,1,0),
        thicknesses=[10], vacuums=[10,20], save_slabs=False, is_symmetric=False, 
        save_metadata=False)

        self.assertEqual(len(sym_true), 2)
        self.assertEqual(len(sym_false), 2)
        self.assertEqual(type(sym_false[0]['slab']), Slab)
    
    def test_save_to_file(self): 
        ytos_slabs = generate_slabs(structure=self.ytos, 
        hkl=(0,0,1), thicknesses=[10], vacuums=[10, 20], save_metadata=False)
        
        self.assertIsNone(ytos_slabs)
        
        # Check the files created 
        self.assertEqual(len(os.listdir('Y2Ti2S2O5')), 1)
        self.assertIn('POSCAR_001_10_10_4.vasp', os.listdir('Y2Ti2S2O5'))

        # Clean up - get rid of directory created 
        shutil.rmtree('Y2Ti2S2O5') 
    
    def test_save_metadata(self):
        ytos_slabs = generate_slabs(structure=self.ytos, 
        hkl=(0,0,1), thicknesses=[10], vacuums=[10, 20], save_slabs=False)
        
        self.assertIsNotNone(ytos_slabs)
        
        # Check the files created 
        self.assertIn('Y2Ti2S2O5_metadata.json', os.listdir(os.getcwd()))

        if os.path.isfile('Y2Ti2S2O5_metadata.json'): 
            os.remove('Y2Ti2S2O5_metadata.json')

    def test_selective_dynamics(self): 
        ytos_slabs = generate_slabs(structure=self.ytos, 
        hkl=(0,0,1), thicknesses=[30,50], vacuums=[20], layers_to_relax=1, 
        save_slabs=False)

        self.assertEqual(ytos_slabs[1]['slab'][28].properties['selective_dynamics'], 
        [0.0, 0.0, 0.0])
        self.assertEqual(ytos_slabs[1]['slab'][0].properties['selective_dynamics'], 
        [1.0, 1.0, 1.0])
        self.assertWarnsRegex(UserWarning, 'Some slabs were too thin to fix the '
        'centre of the slab. Slabs with no selective dynamics applied '
        'are: 001_30_20_4')
        self.assertEqual(type(ytos_slabs[0]['slab']), Slab)
         
        if os.path.isfile('Y2Ti2S2O5_metadata.json'): 
            os.remove('Y2Ti2S2O5_metadata.json')