import unittest

import hotspots.fragment_hotspot_maps as fhm
from ccdc_development.grid import Grid

class Test_BuildLocation(unittest.TestCase):

    def setUp(self):
        self.apolar = Grid.from_file("input_files/apolar.grd")

        indices = (123, 106, 100)
        i = 0
        locations = []
        self.build_location = fhm._BuildLocation( self.apolar, indices, i, mode="volume" )

    def test_apolar_island(self):
        self.assertEquals(self.apolar,self.build_location.apolar_island)

    def test_indicies(self):
        self.assertEquals((123,106,100),self.build_location.indices)

    def test_coordinates(self):
        self.assertAlmostEquals(26.5, self.build_location.coordinates[0],1)


if __name__ == "__main__":
    unittest.main()