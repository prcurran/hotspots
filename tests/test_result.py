import unittest
import numpy as np

from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.result import Results
from hotspots.grid_extension import Grid
from hotspots.protein_extension import Protein
from hotspots.pharmacophore_extension import ProteinPharmacophoreModel
from pprint import pprint


class TestStaticMethods(unittest.TestCase):
    def setUp(self) -> None:
        self.acceptor_grid = Grid.from_file("testdata/6y2g_A/acceptor_mvon_smoothed.grd")

        self.binding_site = Protein.from_file("testdata/6y2g_A/binding_site.pdb")
        self.protein_pharmacophore = ProteinPharmacophoreModel()
        self.protein_pharmacophore.feature_definitions = ["acceptor_projected", "donor_projected", "donor_ch_projected"]
        self.protein_pharmacophore.detect_from_prot(self.binding_site)
        self.donor_features = [f for f in self.protein_pharmacophore.selected_features
                               if f.identifier == "donor_projected"]

    def test_island_to_atom(self):
        self.acceptor_island = self.acceptor_grid.islands(5)[0]
        self.acceptor_island.write("testdata/6y2g_A/first_acceptor_island.grd")

        X = Results._island_to_atom(self.donor_features, self.acceptor_island, tolerance=3)
        print(X)




class TestResults(unittest.TestCase):
    def setUp(self):
        # single result
        self.hotspot_reader = HotspotReader(path="testdata/6y2g_A/out.zip")
        self.result = self.hotspot_reader.read()[0]

        # HotspotReader needs fixing
        self.result.buriedness = Grid.from_file("testdata/6y2g_A/buriedness.grd")

    def test_map_features_to_protein(self):
        pc = self.result._map_features_to_protein()
        print(pc)


if __name__ == '__main__':
    unittest.main()