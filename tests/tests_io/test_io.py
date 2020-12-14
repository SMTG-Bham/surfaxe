import unittest
import warnings
import os
from pathlib import Path
from pymatgen import Structure
from pymatgen.core.surface import Slab
from surfaxe.io import load_config_dict, slab_from_file

class LoadTestCase(unittest.TestCase): 
    def test_load_cd(self): 
        cd1 = load_config_dict('HSE06_config.json')
        cd2 = load_config_dict('waa')
        cd3 = load_config_dict((0,1,2))

        self.assertEqual(cd1['INCAR']['AEXX'], 0.25)
        self.assertEqual(cd1['INCAR']['ALGO'], 'All')
        self.assertEqual(cd2['INCAR']['GGA'], 'PS')
        self.assertEqual(cd2['INCAR']['ALGO'], 'Normal')
        self.assertEqual(cd3['INCAR']['ALGO'], 'Normal')

class SlabFromFileTestCase(unittest.TestCase): 
    def setUp(self): 
        self.slab = str(Path(__file__).parents[2].joinpath('example_data/analysis/POSCAR_LTA_010'))

    def test_slab_from_file(self): 
        slab = slab_from_file(self.slab, (0,1,0))

        self.assertEqual(type(slab), Slab)