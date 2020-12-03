import unittest
import warnings
import os
from pymatgen import Structure
from surfaxe.analysis import simple_nn

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(THIS_DIR, os.pardir, 'data/')

class SimplennTestCase(unittest.TestCase):

    def setUp(self):
        self.poscar = os.path.join(data_dir,'CONTCAR_SnO2')

    def test_simple_nn(self):
        coord_data = simple_nn(start=self.poscar, elements=['Sn','O'], save_csv=False)
        self.assertEqual(coord_data.shape, (90,4))
        self.assertEqual(len(coord_data.loc[coord_data['cn_start'] == 3]), 58)


