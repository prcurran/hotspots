import unittest

import os
from ccdc.protein import Protein
from hotspots.grid_extension import Grid
from hotspots.wrapper_ghecom import Ghecom


class TestGhecom(unittest.TestCase):
    def setUp(self):
        self.protein = Protein.from_file("testdata/1hcl/protein.pdb")
        self.protein.remove_all_waters()
        self.protein.add_hydrogens()
        self.template = Grid.initalise_grid([atm.coordinates for atm in self.protein.atoms])

    def testenvvar(self):
        ghecom_path = f"{os.environ['GHECOM_EXE']}"
        self.assertTrue(os.path.isfile(ghecom_path) and os.access(ghecom_path, os.X_OK))

    def testrun(self):
        self.ghecom = Ghecom()
        self.ghecom.temp = "testdata/wrapper_ghecom"
        self.out = self.ghecom.run(self.protein)
        self.assertTrue(os.path.exists(self.out))

    def testpdb_to_grid(self):
        path = "testdata/wrapper_ghecom/ghecom_out.pdb"
        result = Ghecom.pdb_to_grid(path, self.template)
        result.write("testdata/wrapper_ghecom/test.grd")

        # add some checks


if __name__ == '__main__':
    unittest.main()