import unittest
from hotspots.grid_extension import Grid
from pprint import pprint
from hotspots.wrapper_pymol import PyMOLCommands, PyMOLFile
from ccdc.io import MoleculeReader
import numpy as np


def random_grid(num_of_spheres, return_coords=False, radius=1, value=8, scaling='linear'):

    # something in around the 6Y2G binging site (might be needed later)
    mol = MoleculeReader("testdata/6y2g_A/A_mol.mol2")[0]
    g = Grid.initalise_grid([a.coordinates for a in mol.atoms])

    for i in range(num_of_spheres):
        pnt = [np.random.randint(low=2, high=ax - 2, size=1) for ax in g.nsteps]
        g.set_sphere(point=g.indices_to_point(pnt[0],
                                              pnt[1],
                                              pnt[2]),
                     radius=radius,
                     value=value,
                     scaling=scaling)
    if return_coords:
        return g,
    else:
        return g


class TestGrid(unittest.TestCase):
    def setUp(self):
        # A buriedness grid
        self.buriedness = Grid.from_file("testdata/result/buriedness.grd")
        self.single_peak = random_grid(1)

    def test_neighbourhood(self):
        # check the catchment critera
        neighbours = Grid.neighbourhood(i=0, j=0, k=0, high=self.buriedness.nsteps, catchment=1)
        self.assertEqual(3, len(neighbours))

        neighbours1 = Grid.neighbourhood(i=5, j=5, k=5, high=self.buriedness.nsteps, catchment=1)
        self.assertEqual(6, len(neighbours1))

        f = PyMOLFile()
        for i, n in enumerate(neighbours + neighbours1):
            f.commands += PyMOLCommands.sphere(f"sphere_{i}",
                                               (0, 0, 1, 1),
                                               n,
                                               0.1)
            f.commands += PyMOLCommands.load_cgo(f"sphere_{i}", f"sphere_{i}_obj")
        f.commands += PyMOLCommands.sphere("centre",
                                           (1,0,0, 1),
                                           (5,5,5),
                                           0.1)
        f.commands += PyMOLCommands.load_cgo("centre", "centre_obj")
        f.commands += PyMOLCommands.sphere("centre1",
                                           (1,0,0, 1),
                                           (0,0,0),
                                           0.1)
        f.commands += PyMOLCommands.load_cgo("centre1", "centre1_obj")

        f.write("testdata/grid_extension/neighbourhood.py")

    def test_edge_detection(self):
        self.selected = Grid.from_file("testdata/result/molA.grd")
        edge = self.selected.edge_detection()

        self.assertLess(len(edge), self.selected.count_grid())
        f = PyMOLFile()
        for i, n in enumerate(edge):
            f.commands += PyMOLCommands.sphere(f"sphere_{i}",
                                               (0, 0, 1, 1),
                                               n,
                                               0.1)
            f.commands += PyMOLCommands.load_cgo(f"sphere_{i}", f"sphere_{i}_obj")

        f.write("testdata/grid_extension/edge_detection.py")

    def test_get_peaks(self):
        def visualise(g, gname, peaks):
            # visualise
            gobj_name = f"{gname}_grid"
            surface_objname = f"surface_{gobj_name}"
            vis = f"testdata/grid_extension/{gname}.py"
            gpath = f"testdata/grid_extension/{gname}.grd"

            g.write(gpath)

            f = PyMOLFile()
            f.commands += PyMOLCommands.load(fname=f"{gname}.grd", objname=gobj_name)
            f.commands += PyMOLCommands.isosurface(grd_name=gobj_name, isosurface_name=surface_objname,
                                                   level=5, color="yellow")

            for i, p in enumerate(peaks):
                f.commands += PyMOLCommands.sphere(objname=f"peak_obj_{i}",
                                                   rgba=[1, 1, 1, 1],
                                                   coords=p,
                                                   radius=0.5)
                f.commands += PyMOLCommands.load_cgo(obj=f"peak_obj_{i}", objname=f"peak_{i}")

            f.write(vis)

        # simple case
        self.single_peak = random_grid(1)
        peaks = self.single_peak.get_peaks(min_distance=1)
        self.assertEqual(1, len(peaks))
        # visualise(g=self.single_peak, gname="single_peak", peaks=peaks)

        # single local maxima, multiple pixels
        self.cluster_identical_peak = random_grid(1, scaling='None')
        peaks_1 = self.cluster_identical_peak.get_peaks(min_distance=1)
        self.assertEqual(1, len(peaks_1))
        # visualise(g=self.cluster_identical_peak, gname="cluster_identical", peaks=peaks_1)

        # multiple single pixel local maxima
        self.multi_peak = random_grid(5, scaling='linear')
        peaks_2 = self.multi_peak.get_peaks(min_distance=1)
        self.assertEqual(5, len(peaks_2))
        # visualise(g=self.multi_peak, gname="multi", peaks=peaks_2)

        # multiple multiple pixel local maxima
        self.multi_multi = random_grid(5, scaling='None')
        peaks_3 = self.multi_multi.get_peaks(min_distance=1)
        self.assertEqual(5, len(peaks_3))
        # visualise(g=self.multi_multi, gname="multi_multi", peaks=peaks_2)

if __name__ == '__main__':
    unittest.main()

