import unittest
import warnings
import os
from pathlib import Path
from pymatgen import Structure
from surfaxe.analysis import simple_nn, complex_nn

data_dir = str(Path(__file__).parents[1].joinpath('example_data/analysis'))

class NNTestCase(unittest.TestCase):

    def setUp(self):
        self.SnO2 = os.path.join(data_dir, 'CONTCAR_SnO2')
        self.lta = os.path.join(data_dir,'POSCAR_LTA_010')

    def test_simple_nn(self):
        coord_data = simple_nn(start=self.SnO2, elements=['Sn','O'], save_csv=False)
        self.assertEqual(coord_data.shape, (90,4))
        self.assertEqual(len(coord_data.loc[coord_data['cn_start'] == 3]), 58)

    def test_complex_nn(self):
        coord_data = complex_nn(start=self.lta, 
           elements=['La', 'Ag', 'Ti', 'O', 'S'], 
           cut_off_dict={('Ag','S'): 3.09, ('La','O'): 2.91,
                          ('La','S'): 3.559, ('Ti','O'): 2.35,
                          ('Ti','S'): 2.91,}, save_csv=False)
        self.assertEqual(coord_data.shape, (240,4))
        



