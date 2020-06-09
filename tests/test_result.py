import unittest
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.result import Results
from pprint import pprint


class TestResults(unittest.TestCase):
    def setUp(self):
        # single result
        self.hotspot_reader = HotspotReader(path="testdata/result/out.zip")
        self.result = self.hotspot_reader.read()

    def test_map_features_to_protein(self):
        self.result._map_features_to_protein()
