import unittest
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.result import Results
from pprint import pprint


class TestHotspotReader(unittest.TestCase):
    def setUp(self):
        # single result
        self.hotspot_reader = HotspotReader(path="testdata/hs_io/out.zip")

    def test_read(self):
        self.result = self.hotspot_reader.read()
        assert isinstance(self.result, Results)


class TestHotspotWriter(unittest.TestCase):
    def setUp(self):
        # single result
        self.result = HotspotReader(path="testdata/hs_io/out.zip").read()

    def test_get_labels(self):
        labs = self.result.grid_labels()

        pprint(labs)

    def test_write(self):
        with HotspotWriter("testdata/hs_io/test_write") as w:
            w.write(self.result)



