import unittest
import numpy as np
from hotspots.hs_io import HotspotReader
from hotspots.grid_extension import Grid
from hotspots.protein_extension import Protein, centroid
from hotspots.wrapper_pymol import PyMOLCommands, PyMOLFile
from ccdc import io

class TestBindingSiteFromGrid(unittest.TestCase):
    def setUp(self):
        self.prot = Protein.from_file("testdata/result/binding_site.pdb")
        self.grid = Grid.from_file("testdata/result/molA.grd")

    def test_centroid(self):
        cent = centroid(np.array([a.coordinates for a in self.prot.residues[0].atoms]))

        expected = (1.889772727272726, 2.9687272727272713, -21.237954545454546)
        self.assertAlmostEqual(cent[0], expected[0])
        self.assertAlmostEqual(cent[1], expected[1])
        self.assertAlmostEqual(cent[2], expected[2])

    def test_detect_from_grid(self):
        bs = Protein.BindingSiteFromGrid._detect_from_grid(self.prot, self.grid, 4)
        self.assertEqual(18, len(bs))

        # for visual inspection
        f = PyMOLFile()
        f.commands += PyMOLCommands.load("binding_site.pdb", "abc")

        for res in bs:
            f.commands += PyMOLCommands.select("sele", f'resi {res.identifier.split(":")[1][3:]}')
            f.commands += PyMOLCommands.show("sticks", "sele")

        f.write("testdata/protein_extension/test_bindingsitefromgrid.py")

    def test_BindingSiteFromGrid(self):

        bs = Protein.BindingSiteFromGrid(self.prot, self.grid, within=6)
        print(len(bs.residues))
        self.assertEqual(28, len(bs.residues))

        f = PyMOLFile()
        f.commands += PyMOLCommands.load("binding_site.pdb", "abc")

        with io.MoleculeWriter("testdata/protein_extension/binding_site.pdb") as w:
            w.write(bs.protein)

        for res in bs.residues:
            f.commands += PyMOLCommands.select("sele", f'resi {res.identifier.split(":")[1][3:]}')
            f.commands += PyMOLCommands.show("sticks", "sele")

        f.write("testdata/protein_extension/test_bindingsitefromgrid.py")




