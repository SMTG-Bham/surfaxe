import unittest
from pathlib import Path
from pymatgen.core.surface import Slab
from surfaxe.io import _load_config_dict, slab_from_file

class LoadTestCase(unittest.TestCase): 
    def test_load_cd(self): 
        cd1 = _load_config_dict('HSE')
        cd2 = _load_config_dict('waa')
        cd3 = _load_config_dict((0,1,2))

        self.assertEqual(cd1['INCAR']['AEXX'], 0.25)
        self.assertEqual(cd1['INCAR']['ALGO'], 'All')
        self.assertEqual(cd2, cd3)
        self.assertNotEqual(cd1, cd2)
        self.assertEqual(cd2['INCAR']['GGA'], 'PS')
        self.assertEqual(cd3['INCAR']['ALGO'], 'Normal')
        self.assertRaises(FileNotFoundError, _load_config_dict, 'HSE06.json')

class SlabFromFileTestCase(unittest.TestCase): 
    def setUp(self): 
        self.slab = str(Path(__file__).parents[2].joinpath('example_data/analysis/POSCAR_LTA_010'))

    def test_slab_from_file(self): 
        slab = slab_from_file(self.slab, (0,1,0))

        self.assertEqual(type(slab), Slab)