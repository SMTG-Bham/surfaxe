import unittest
import warnings
import os
from pathlib import Path
from pymatgen import Structure
from surfaxe.analysis import simple_nn, complex_nn, cart_displacements, \
bond_analysis, electrostatic_potential

data_dir = str(Path(__file__).parents[2].joinpath('example_data/analysis'))

class NNTestCase(unittest.TestCase):

    def setUp(self):
        self.SnO2 = os.path.join(data_dir, 'CONTCAR_SnO2')
        self.lta = os.path.join(data_dir,'POSCAR_LTA_010')
        self.lta_end = os.path.join(data_dir, 'CONTCAR_LTA_010')

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
        self.assertEqual(coord_data['nn_start'][52], 'O O S S S')

        end_data = complex_nn(start=self.lta, 
           elements=['La', 'Ag', 'Ti', 'O', 'S'], 
           cut_off_dict={('Ag','S'): 3.09, ('La','O'): 2.91,
                          ('La','S'): 3.559, ('Ti','O'): 2.35,
                          ('Ti','S'): 2.91,}, 
           save_csv=False, end=self.lta_end)
        self.assertEqual(end_data.shape, (240, 6))
        self.assertNotEqual(end_data['nn_start'][101], end_data['nn_end'][101])
        

class CartDisplacementsTestCase(unittest.TestCase): 

    def setUp(self): 
        self.start = os.path.join(data_dir, 'POSCAR_LTA_010')
        self.end = os.path.join(data_dir, 'CONTCAR_LTA_010')
    
    def test_cart_displacements(self): 
        cart_data = cart_displacements(start=self.start, 
            end=self.end, elements=['La', 'Ti', 'Ag', 'S', 'O'],
            save_txt=False)
        self.assertIsNotNone(cart_data)
        self.assertEqual(len(cart_data['site']), 192)

class BondAnalysisTestCase(unittest.TestCase): 

    def setUp(self): 
        self.structure = os.path.join(data_dir, 'CONTCAR_SnO2')

    def test_bond_analysis(self): 
        bonds_data = bond_analysis(structure=self.structure, 
        bonds = [['Sn', 'O']], save_csv=False, save_plt=False)
        self.assertEqual(bonds_data.shape, (30,3))

class ElectrostaticPotentialTestCase(unittest.TestCase): 

    def setUp(self): 
        self.locpot = os.path.join(data_dir, 'LOCPOT')

    def test_electrostatic_potential(self): 
        potential_data = electrostatic_potential(10, locpot=self.locpot, 
        save_csv=False, save_plt=False)
        self.assertEqual(potential_data.shape, (1372,2))
        self.assertEqual(potential_data['planar'][394], -10.668138966414821)