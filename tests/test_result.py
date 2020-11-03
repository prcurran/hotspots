import unittest
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.result import Extractor
from hotspots.grid_extension import Grid
from ccdc.io import MoleculeReader
from ccdc.utilities import PushDir
import os


class TestExtractor(unittest.TestCase):
    def setUp(self) -> None:
        with HotspotReader(path="testdata/result/data/out.zip") as r:
            self.result = r.read()

        for p, g in self.result.super_grids.items():
            self.result.super_grids[p] = g.dilate_by_atom()
        # bin
        self.bin = "testdata/result/Extractor/bin"

        # reuse
        self.out = "testdata/result/Extractor"

    def testconstruction(self):
        extractor = Extractor(self.result)

        # extractor.single_grid.write(os.path.join(self.out, "2vta_single_grid.grd"))

        hr = extractor.extract_volume()

        with HotspotWriter(self.bin) as w:
            w.write(hr)


class TestResult(unittest.TestCase):
    def testscore_atoms_as_spheres(self):
        with PushDir("testdata/result/data"):
            mols = [m for m in MoleculeReader("gold_docking_poses.sdf")]

            # create a grid which can contain all docking poses
            small_blank = Grid.initalise_grid(coords={atm.coordinates for mol in mols for atm in mol.heavy_atoms},
                                              padding=2)

            # read hotspot maps
            with HotspotReader(path="out.zip") as r:
                self.result = r.read()

            # dilate the grids
            for p, g in self.result.super_grids.items():
                self.result.super_grids[p] = g.dilate_by_atom()

            # shrink hotspot maps to save time
            sub_grids = {p: Grid.shrink(small=small_blank, big=g)
                         for p, g in self.result.super_grids.items()}

            # create single grid
            mask_dic, sg = Grid.get_single_grid(sub_grids)

            self.result.super_grids = mask_dic

            # set background to 1
            self.result.set_background()
            self.result.normalize_to_max()

            print([g.extrema for p, g in self.result.super_grids.items()])

            for m in mols[:1]:
                s = self.result.score_atoms_as_spheres(m, small_blank)
                print(s)

    def test_docking_constraint_atoms(self):
        with PushDir("testdata/result/data"):
            # read hotspot maps
            with HotspotReader(path="out.zip") as r:
                self.result = r.read()

            print(self.result._docking_constraint_atoms())



if __name__ == '__main__':
    unittest.main()