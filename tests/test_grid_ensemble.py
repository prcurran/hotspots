from __future__ import print_function, division

import unittest
from hotspots.grid_extension import Grid, _GridEnsemble
import numpy as np
from glob import glob
import os
from os.path import join, exists
import shutil


class TestEnsembleSyntheticData(unittest.TestCase):

    @staticmethod
    def make_test_data():
        def make_grid(offset, vals, idxs, nsteps):
            grid_origin = offset
            grid_far_corner = (
            offset[0] + (nsteps[0] - 1) * 0.5, offset[1] + (nsteps[1] - 1) * 0.5, offset[2] + (nsteps[2] - 1) * 0.5)
            out_grid = Grid(origin=grid_origin,
                            far_corner=grid_far_corner,
                            spacing=0.5,
                            _grid=None,
                            default=0)
            for (nx, ny, nz), v in zip(idxs, vals):
                # print(nx, ny, nz, v)
                # print(int(nx-offset[0]*2), int(ny-offset[1]*2), int(nz-offset[2]*2))
                out_grid.set_value(int(nx - offset[0] * 2), int(ny - offset[1] * 2), int(nz - offset[2] * 2), v)
            return out_grid

        # Set up the numpy array:
        nsteps = (40, 50, 45)

        # Fill with a thousand random floats between 0 and 40 (roughly hotspots range)
        np.random.seed(3)
        values = np.random.uniform(1, 40, 1000)
        ind_x = np.random.randint(5, nsteps[0] - 5, size=1000)
        ind_y = np.random.randint(5, nsteps[1] - 5, size=1000)
        ind_z = np.random.randint(5, nsteps[2] - 5, size=1000)

        # Create 10 grids, with an origin offset of up to 10 grid points in each direction
        grid_spacing = 0.5
        offset_x = np.random.randint(-5, 5, 10) * grid_spacing
        offset_y = np.random.randint(-5, 5, 10) * grid_spacing
        offset_z = np.random.randint(-5, 5, 10) * grid_spacing

        offsets = zip(offset_x, offset_y, offset_z)
        print(offsets)
        indices = zip(ind_x, ind_y, ind_z)

        tmp_dir = join(os.getcwd(), "tmp_gridensemble_test")
        if not exists(tmp_dir):
            os.mkdir(tmp_dir)

        grid_list = []
        for i in range(len(offsets)):
            off = offsets[i]
            g = make_grid(off, values, indices, nsteps)
            print(g.nsteps)
            grid_list.append(g)
            g.write(join(tmp_dir, "test_apolar_{}.ccp4".format(str(i))))

        return grid_list

    def setUp(self):
        self.grid_list = self.make_test_data()
        self.tmp_dir =  join(os.getcwd(), "tmp_gridensemble_test")
        self.grid_paths = glob(join(self.tmp_dir, "test_apolar_*.ccp4"))
        print(self.grid_paths)
        self.grid_ensemble = _GridEnsemble()
        self.grid_ensemble = _GridEnsemble()
        self.max_grid = self.grid_ensemble.from_grid_list(self.grid_list, os.getcwd(), "test", "apolar")
        self.values = self.grid_ensemble.results_array[self.grid_ensemble.results_array.nonzero()]

    def test_from_hotspot_maps(self):
        # hard coded because 988 values in test set - perhaps there's a better way to test this?
        grid_ensemble = _GridEnsemble()
        max_grid = grid_ensemble.from_hotspot_maps(self.grid_paths, os.getcwd(), "test", "apolar")
        values = grid_ensemble.results_array[grid_ensemble.nonzeros]
        self.assertEqual(len(values),
                         988,
                         msg="Check number of values")
        self.assertEqual(grid_ensemble.tup_max_length,
                         len(self.grid_paths),
                         msg="Check tup_max_length assigned")

    def test_from_grid_list(self):
        grid_ensemble = _GridEnsemble()
        max_grid = grid_ensemble.from_grid_list(self.grid_list, os.getcwd(), "test", "apolar")
        values = grid_ensemble.results_array[grid_ensemble.results_array.nonzero()]
        self.assertEqual(len(values),
                         988,
                         msg="Check number of values")
        self.assertEqual(grid_ensemble.tup_max_length,
                         len(self.grid_list),
                         msg="Check tup_max_length assigned")

    def test_get_gridpoint_max(self):

        maxes = self.grid_ensemble.get_gridpoint_max()
        self.assertEqual(len(maxes),len(self.values), msg="Check get_max has correct number of values")
        arr_max = np.array(maxes)
        arr_vals = np.array([self.values[i][0] for i in range(len(self.values))])
        self.assertTrue(np.array_equal(arr_max, arr_vals), msg="Check get_max")

    def test_get_gridpoint_means(self):
        means = self.grid_ensemble.get_gridpoint_means()
        self.assertEqual(len(means), len(self.values), msg="Check get_max has correct number of values")
        arr_mean = np.array(means)
        arr_vals = np.array([self.values[i][0] for i in range(len(self.values))])
        self.assertTrue(np.allclose(arr_mean, arr_vals), msg="Check get_means")

    def test_get_gridpoint_means_spread(self):
        means = self.grid_ensemble.get_gridpoint_means_spread()
        self.assertIsInstance(means, list)
        self.assertEqual(len(means), len(self.values)*10)

    def test_get_gridpoint_ranges(self):
        ranges = self.grid_ensemble.get_gridpoint_ranges()
        self.assertIsInstance(ranges, list)
        self.assertEqual(len(set(ranges)), 1)

    def test_output_grid(self):
        #g_0 = Grid.from_file(self.grid_paths[0])
        g_0 = self.grid_list[0]
        print(g_0.count_grid())
        print(self.max_grid.count_grid())
        max_g, ref_g = Grid.common_grid([self.max_grid, g_0])
        diff_grid = (max_g - ref_g)
        nx, ny, nz = diff_grid.nsteps
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    if diff_grid.value(x,y,z)!=0:
                        print(diff_grid.value(x,y,z))
        self.assertEqual((max_g-ref_g).count_grid(), 0, msg="Testing the max_grid")

        means_grid = self.grid_ensemble.output_grid(mode="mean", save=False)
        mean_g, ref_g = Grid.common_grid([means_grid, g_0])
        self.assertEqual((max_g - ref_g).count_grid(), 0, msg="Testing the means_grid")

        ranges_grid = self.grid_ensemble.output_grid(mode="ranges", save=False)
        self.assertEqual(ranges_grid.count_grid(), 0, msg="Testing the ranges grid")

        other_g = self.grid_ensemble.output_grid(mode="bla", save=False)
        self.assertIsNone(other_g)

    def tearDown(self):
        if exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        #os.remove("test_max_apolar.ccp4")



if __name__ == "__main__":
    unittest.main()