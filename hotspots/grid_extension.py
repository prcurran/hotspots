'''
The main class of the :mod:`hotspots.grid_extension.Grid`.

This is an internal extension of :class:`ccdc.grid.Grid` that adds in potential new _features
for review and development. For example, gaussian smoothing function:

.. code-block:: python
#
# >>> from ccdc_development import Grid
# >>> from scipy import ndimage
# >>> import numpy as np
# >>>
# >>> grd = Grid.from_file(<path/to/file>)
# >>> new_grd = grd.gaussian(sigma=0.3)
# >>>

I.e. this returns a new grid object in which gaussian smoothing has been applied
'''
#############################################################################
from __future__ import print_function, division

import collections
import operator

import numpy as np
from ccdc import utilities
from hotspots.hs_utilities import Helper
from scipy import ndimage
from skimage import feature
from os.path import join, basename
from scipy.stats import norm
import matplotlib.pyplot as plt
import pickle


Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


class Grid(utilities.Grid):
    """
    Initialisation handled in the main class
    """
    def grid_values(self, threshold=0):
        """"""
        array = self.get_array()
        masked_array = np.ma.masked_less_equal(array, threshold)
        return masked_array.compressed()

    def grid_score(self, threshold=0, percentile=75):
        """
        take a group and return average score of points above threshold
        :param g:
        :param threshold
        :param percentile

        :return:
        """
        array = self.get_array()
        masked_array = np.ma.masked_less_equal(array, threshold)
        values = masked_array.compressed()

        if len(values) == 0:
            return 0
        else:
            return np.percentile(values, percentile)

    def indices_to_point(self, i, j, k):
        """
        Return x,y,z coordinate for a given grid index

        :param i: int, indice (x-axis)
        :param j: int, indice (y-axis)
        :param k: int, indice (z-axis)
        :param g: a :class: `ccdc.utilities.Grid` instance
        :return: float(x), float(y), float(z)
        """

        ox, oy, oz, = self.bounding_box[0]
        gs = self.spacing
        return ox + float(i) * gs, oy + float(j) * gs, oz + gs * float(k)

    def point_to_indices(self, p):
        """
        Return the nearest grid index for a given point

        :param p: tup, (float(x), float(y), float(z)), a coordinate on a 3D grid
        :param g: a :class: `ccdc.utilities.Grid` instance
        :return: int(x), int(y), int(z)
        """

        gs = self.spacing
        rx, ry, rz = [round(i / gs) for i in p]
        ox, oy, oz = [round(i / gs) for i in self.bounding_box[0]]
        return int(rx - ox), int(ry - oy), int(rz - oz)

    def gaussian(self, sigma=0.2):
        """
        gaussian smoothing function, method of reducing noise in output
        :param sigma:
        :return:
        """
        s = (sigma, sigma, sigma, 0)
        nx, ny, nz = self.nsteps
        scores = np.zeros((nx, ny, nz, 1))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    scores[i, j, k, 0] += self.value(i, j, k)
        smoothed = ndimage.filters.gaussian_filter(scores, sigma=s)
        grid = self.copy_and_clear()
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    grid.set_value(i, j, k, smoothed[i, j, k, 0])
        return grid

    def contains_point(self, point, threshold=0, tolerance=0):
        """
        determines whether a set of coordinates are within a grids boundary
        :param point:
        :return:
        """
        mini = self.bounding_box[0]
        maxi = self.bounding_box[1]
        if self.value_at_point(point) >= threshold:
            return all([mini.x - tolerance < point[0] < maxi.x + tolerance,
                        mini.y - tolerance < point[1] < maxi.y + tolerance,
                        mini.z - tolerance < point[2] < maxi.z + tolerance])
        else:
            return False

    def get_array(self):
        """
        convert grid object to np.array
        :param self:
        :return:
        """
        nx, ny, nz = self.nsteps
        array = np.zeros((nx, ny, nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    array[i, j, k] += self.value(i, j, k)
        return array

    def restricted_volume(self, volume=75):
        """
        returns a grid with of a defined volume
        :param self:
        :param volume:
        :return:
        """
        grid = self.copy_and_clear()
        max_points = int(float(volume) / 0.125)
        nx, ny, nz = self.nsteps
        rank_dict = {}

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    value = float(self.value(i, j, k))
                    if value in rank_dict:
                        rank_dict[value].append((i, j, k))
                    else:
                        rank_dict.update({value: [(i, j, k)]})

        top_points = sorted((float(x) for x, y in rank_dict.iteritems()), reverse=True)
        indices = [pts for key in top_points for pts in rank_dict[key]]

        for i in indices[:max_points]:
            grid.set_value(i[0], i[1], i[2], self.value(i[0], i[1], i[2]))

        return Grid.super_grid(1, *grid.islands(threshold=1))

    def centroid(self):
        """
        returns centre of grid
        :param self:
        :return:
        """
        return ((self.bounding_box[0][0] + self.bounding_box[1][0])/2,
                (self.bounding_box[0][1] + self.bounding_box[1][1])/2,
                (self.bounding_box[0][2] + self.bounding_box[1][2])/2
                )

    def deduplicate(self, major, threshold=12, tolerance=2):
        """
        method to deduplicate two grids
        :param self: clearing grid
        :param major: overriding grid
        :param threshold: island threshold
        :param tolerance: tolerance in contains point metho
        :return:
        """
        if self.bounding_box[0] != major.bounding_box[0] or self.bounding_box[1] != major.bounding_box[1]:
            self = major.common_boundaries(self)

        all_islands = set([jsland for jsland in self.islands(threshold=threshold)])
        bin_islands = set([jsland for jsland in all_islands
                           for island in major.islands(threshold=threshold)
                           if jsland.contains_point(island.centroid(), tolerance=tolerance)
                           or jsland.count_grid() <= 8
                           or Helper.get_distance(jsland.centroid(), island.centroid()) < 4])

        retained_jslands = list(all_islands - bin_islands)

        if len(retained_jslands) == 0:
            blank = major.copy_and_clear()
            return blank
        else:
            temp = Grid.super_grid(0, *retained_jslands)
            blank = self.copy_and_clear()
            return blank.common_boundaries(temp)

    def copy_and_clear(self):
        """
        make a new empty grid
        :param self:
        :return:
        """
        g = self.copy()
        g *= 0
        return g

    def common_boundaries(self, grid):
        """
        Expands parsed grid to the size of self grid (supplied grid should be smaller than self)
        :param grid:
        :return:
        """
        reference = Grid.super_grid(0, self, grid)
        blank = reference.copy_and_clear()

        return Grid.super_grid(0, blank, grid)

    def multi_max_mask(self, grids):
        """
        given 1 or more grid, the value selected is the grid point across the input grids which contains the highest val
        this doesnt make sense
        :param grid:
        :return:
        """
        max_grids = [self > g for g in grids]
        blank = -self.copy_and_clear()
        return reduce(operator.__and__, max_grids, blank)

    def get_best_island(self, threshold, mode="count", peak=None):
        """
        returns the best grid island. Mode: "count" or "score"
        :param threshold: island threshold
        :param mode:
                    -"count" : returns island with most grid points above threshold
                    -"score" : returns island with the largest sum of all grid points over threshold

        :return:
        """
        islands = self.islands(threshold)
        if len(islands) == 0:
            return None

        else:
            island_by_rank = {}
            if mode == "count":
                for island in islands:
                    if peak:
                        if island.contains_point(peak):
                            g = (island > threshold) * island
                            rank = g.count_grid()
                            island_by_rank.update({rank: island})
                        else:
                            continue
                    else:
                        g = (island > threshold) * island
                        rank = g.count_grid()
                        island_by_rank.update({rank: island})

            elif mode == "score":
                for island in islands:
                    if peak:
                        if island.contains_point(peak):
                            nx, ny, nz = island.nsteps
                            island_points = [island.value(i, j, k)
                                             for i in range(nx) for j in range(ny) for k in range(nz)
                                             if island.value(i, j, k) >= threshold]
                            rank = sum(island_points)
                            island_by_rank.update({rank: island})
                        else:
                            continue
                    else:
                        nx, ny, nz = island.nsteps
                        island_points = [island.value(i, j, k)
                                         for i in range(nx) for j in range(ny) for k in range(nz)
                                         if island.value(i, j, k) >= threshold]
                        rank = sum(island_points)
                        island_by_rank.update({rank: island})

            else:
                raise IOError("mode not supported")

            if len(island_by_rank) == 0:
                return None
            else:
                rank = sorted(island_by_rank.keys(), reverse=True)[0]
                return island_by_rank[rank]

    def minimal(self):
        """Reduces grid size to the minimal dimensions"""
        return Grid.super_grid(1, *self.islands(threshold=1))

    def limit_island_size(self, npoints, threshold=10):
        """for a given grid, there are no islands above npoints (at any value)"""
        g = (self > 10) * self
        all_islands = []
        for island in g.islands(threshold):
            if island.count_grid() > npoints:
                all_islands.append(island.top_points(npoints=npoints))
            else:
                all_islands.append(island)
        return Grid.super_grid(0, *all_islands)

    def top_points(self, npoints):
        """orders points and returns the top npoints"""
        pts = {}
        nx, ny, nz = self.nsteps

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    val = self.value(i, j, k)

                    if val in pts:
                        pts[val].append((i, j, k))
                    else:
                        pts.update({self.value(i, j, k): [(i, j, k)]})

        sorted_pts = sorted(pts.items(), key=lambda x: x[0], reverse=True)

        thres = self._get_threshold(sorted_pts, npoints=npoints)
        #thres = sorted_pts[npoints][0]

        return (self > thres) * self

    def step_out_mask(self, nsteps=2):
        """Add one step in all directions to Grid boundary, returns blank grid"""
        origin = (self.bounding_box[0].x - (self.spacing * nsteps),
                  self.bounding_box[0].y - (self.spacing * nsteps),
                  self.bounding_box[0].z - (self.spacing * nsteps))

        far_corner = (self.bounding_box[1].x + (self.spacing * nsteps),
                      self.bounding_box[1].y + (self.spacing * nsteps),
                      self.bounding_box[1].z + (self.spacing * nsteps))

        return Grid(origin=origin,
                    far_corner=far_corner,
                    spacing=0.5,
                    default=0,
                    _grid=None)

    @staticmethod
    def _get_threshold(sorted_points, npoints):
        """private method, hands off"""
        count = []
        for value, pts in sorted_points:
            count.extend(pts)
            if len(count) >= npoints:
                return value
            else:
                continue

    @staticmethod
    def from_array(fname):
        """
        creates a grid from array
        :param fname:
        :return:
        """
        grid = Grid(origin=[-35.00, -42.00, 44.00],
                    far_corner=[59.00, 53.00, 54.00],
                    spacing=0.5,
                    default=0.0,
                    _grid=None)

        array = np.load(fname)
        indices = np.nonzero(array)
        values = array[indices]
        as_triads = zip(*indices)

        for (i, j, k), v in zip(as_triads, values):
            grid._grid.set_value(int(i), int(j), int(k), v)

        return grid

    @staticmethod
    def common_grid(grid_list, padding=1):
        """
        returns two grid with common boundaries
        :param grid_a:
        :param grid_b:
        :param padding:
        :return:
        """
        sg = Grid.super_grid(padding, *grid_list)
        out_g = sg.copy()
        out_g *= 0
        out = [Grid.super_grid(padding, g, out_g) for g in grid_list]
        return out

    def inverse_single_grid(self, mask_dic):
        """
        from a single grid, uses mask_dic to return separated grids divide by interaction type
        :param mask_dic: dict
        :return: dict
        """
        sg = mask_dic["apolar"].common_boundaries(self)
        return {probe: grid * (grid & sg) for probe, grid in mask_dic.items()}

    @staticmethod
    def get_single_grid(grd_dict, mask=True):
        """
        Combines a dictionary of identifier (str): Grid (ccdc.utilties.Grid) to a single grid.
        Grid points of the single_grid are set to the maximum value at each point across all the input grids
        :param grd_dict:
        :param mask:
        :return:
        """
        mask_dic = {}

        if len(set([g.bounding_box for g in grd_dict.values()])) > 1:
            grd_dict = dict(zip(grd_dict.keys(), Grid.common_grid(grd_dict.values())))

        for probe, grid in grd_dict.items():
            other_grids = [grd_dict[p] for p in grd_dict.keys() if p != probe]
            mask_dic.update({probe: grid * grid.multi_max_mask(other_grids)})

        o = [g.bounding_box[0] for g in grd_dict.values()]
        origin = [min([coord.x for coord in o]),
                  min([coord.y for coord in o]),
                  min([coord.z for coord in o])
                  ]

        f = [g.bounding_box[1] for g in grd_dict.values()]
        far_corner = [max([coord.x for coord in f]),
                      max([coord.y for coord in f]),
                      max([coord.z for coord in f])
                      ]

        blank = Grid(origin=origin, far_corner=far_corner, spacing=0.5, default=0, _grid=None)

        if mask:
            return mask_dic, reduce(operator.add, mask_dic.values(), blank)
        else:
            return reduce(operator.add, mask_dic.values(), blank)

    def percentage_overlap(self, other):
        """
        find the percentage overlap of this grid with other.
        :param other: `hotspots.grid_extension.Grid`
        :return:`hotspots.grid_extension.Grid`
        """
        g, h = Grid.common_grid(grid_list=[self, other], padding=1)
        vol = (g > 0).count_grid()
        overlap = (g & h).count_grid()
        return (overlap / vol) * 100

    @staticmethod
    def from_molecule(mol, scaling=1):
        """
        generate a molecule mask where gp within the vdw radius of the molecule heavy atoms are set to 1.0
        :param mol: `ccdc.molecule.Molecule`
        :param padding: int
        :param scaling: float
        :return: `hotspots.grid_extension.Grid`
        """
        g = Grid.initalise_grid(coords=mol.atoms)
        for a in mol.heavy_atoms:
            g.set_sphere(point=a.coordinates,
                         radius=a.vdw_radius * scaling,
                         value=1,
                         scaling='None')
        return g > 0.1

    @staticmethod
    def initalise_grid(coords, padding=1):
        """
        creates a fresh grid using a list of coordinates to define the grid limits
        :param coords: list
        :param padding: int, extra spacing added to grid limits
        :return:
        """
        x = set()
        y = set()
        z = set()
        for c in coords:
            x.add(c.x)
            y.add(c.y)
            z.add(c.z)

        origin = Coordinates(x=round(min(x) - padding),
                             y=round(min(y) - padding),
                             z=round(min(z) - padding))

        far_corner = Coordinates(x=round(max(x) + padding),
                                 y=round(max(y) + padding),
                                 z=round(max(z) + padding))

        return Grid(origin=origin, far_corner=far_corner, spacing=0.5, default=0, _grid=None)

    @staticmethod
    def grow(inner, template, percentile=60):
        """
        experimental
        Dilates grid to the points in the top percentile of the template
        :param template:
        :return:
        """
        expand = inner.max_value_of_neighbours() > 0.1   # remove very small values
        outer = expand.__sub__(inner) * template
        threshold = np.percentile(a=outer.grid_values(threshold=1), q=int(percentile))
        return inner.__add__(outer > threshold)

    def get_peaks(self, min_distance=6, cutoff=2):
        """
        find peak coordinates in grid
        :return:
        """
        peaks = feature.peak_local_max(self.get_array(),
                                       min_distance=min_distance,
                                       threshold_abs=cutoff)
        peak_by_value = {}
        for peak in peaks:
            val = self.value(int(peak[0]), int(peak[1]), int(peak[2]))
            if val > cutoff:
                if val in peak_by_value:
                    peak_by_value[val].append((peak[0], peak[1], peak[2]))
                else:
                    peak_by_value.update({val: [(peak[0], peak[1], peak[2])]})

        average_peaks = []
        for key in peak_by_value.keys():
            x = [point[0] for point in peak_by_value[key]]
            y = [point[1] for point in peak_by_value[key]]
            z = [point[2] for point in peak_by_value[key]]
            average_peaks.append(self.indices_to_point(int(sum(x) / len(x)),
                                                       int(sum(y) / len(y)),
                                                       int(sum(z) / len(z))
                                                       )
                                 )
        return average_peaks

    def value_at_coordinate(self, coordinates, tolerance=1, position=True):
        """
        Uses Grid.value() rather than Grid.value_at_point(). Chris Radoux reported speed issues.
        :param coordinates:
        :param tolerance:
        :return:
        """
        i, j, k = self.point_to_indices(coordinates)
        nx, ny, nz = self.nsteps
        scores = {}

        for di in range(-tolerance, +tolerance + 1):
            for dj in range(-tolerance, +tolerance + 1):
                for dk in range(-tolerance, +tolerance + 1):
                    if 0 < (i + di) < nx and 0 < (j + dj) < ny and 0 < (k + dk) < nz:
                        scores.update({self.value(i + di, j + dj, k + dk): (i + di, j + dj, k + dk)})

        if len(scores) > 0:
            score = sorted(scores.keys(), reverse=True)[0]

            if score < 0.1:
                score = 0
                point = (0, 0, 0)

            else:
                a, b, c = scores[score]
                point = self.indices_to_point(a, b, c)

        else:
            score = 0
            point = (0, 0, 0)

        if position:
            return score, point

        else:
            return score

utilities.Grid = Grid


class _GridEnsemble(object):
    """Class that handles a numpy array of tuples from Hotspot maps of a given probe type across multiple structures"""

    def __init__(self):
        self.prot_name = None
        self.probe = None
        self.path_list = None
        self.tup_max_length = None
        self.results_array = None
        self.common_grid_origin = None
        self.common_grid_far_corner = None
        self.common_grid_nsteps = None
        self.nonzeros = None
        self.out_dir = None
        self.spacing = None

    ###### Functions that generate the ensemble and set attributes

    def _common_grids_from_grid_list(self, grid_list, fnames=None):
        """
        Gets the coordinates of a grid that can fit all other grids (to use as basis for self.results_array)
        :param grid_list: list of 'hotspots.grid_extension.Grid' objects
        :return: list of 'hotspots.grid_extension.Grid' objects
        """
        print("Making array grid {} {}".format(self.prot_name, self.probe))
        common_grids = Grid.common_grid(grid_list)
        
        
        assert (len(set([c.bounding_box for c in common_grids])) == 1), "Common grids don't have same frames"
        self.spacing = grid_list[0].spacing
        # assert (common_grids[0].spacing == 0.5), "Grid spacing not 0.5"
        self.tup_max_length = len(grid_list)
        self.common_grid_origin = common_grids[0].bounding_box[0]
        self.common_grid_far_corner = common_grids[0].bounding_box[1]
        self.common_grid_nsteps = common_grids[0].nsteps

        if fnames:
            self.path_list = fnames

        return common_grids

    def _common_grids_from_paths(self):
        """
        Gets the coordinates of a grid that can fit all other grids (to use as basis for self.results_array)
        :return: list of 'hotspots.grid_extension.Grid' objects
        """
        print('Making array grid')
        grid_list = [Grid.from_file(f) for f in self.path_list if self.probe in basename(f)]
        common_grids = self._common_grids_from_grid_list(grid_list)

        return common_grids

    def _get_results_array(self, common_grids):
        """
        constructs the numpy array containing the list of tuples.
        Adjusts the coordinates based on the origins of the grids
        self.results_array = numpy array
        :return: 
        """
        print('Making results array')

        results_array = np.zeros(self.common_grid_nsteps, dtype=tuple)
        # rec_spacing = 1 / self.spacing

        nx, ny, nz = self.common_grid_nsteps
        for c in common_grids:
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        if c.value(x, y, z) != 0:
                            if isinstance(results_array[x][y][z], tuple):
                                results_array[x][y][z] += (c.value(x, y, z),)
                            else:
                                results_array[x][y][z] = (c.value(x, y, z),)

        self.results_array = results_array
        self.nonzeros = self.results_array.nonzero()
        
    def get_alternative_results_array(self,grid_list):
        """"
        Return to old way of making the array due to bad allocs with Grid.common_grid() method
        """
        # Find bounding box of smallest grid that fits all the grids
        dimensions = np.array([g.bounding_box for g in grid_list])

        self.common_grid_origin = tuple(np.min(dimensions[:, 0, :], axis=0))
        self.common_grid_far_corner = tuple(np.max(dimensions[:, 1, :], axis=0))
        
        origins = dimensions[:, 0, :]
        far_corners = dimensions[:, 1, :]
        if self.spacing:
            rec_spacing = 1 / self.spacing
        else:
            rec_spacing = 2
            
        origin_diff = (origins - self.common_grid_origin)*rec_spacing
        
        if (origins == origins[0]).all() and (far_corners == far_corners[0]).all():
            print("Grids have same dimensions")
            self.common_grid_nsteps = grid_list[0].nsteps
            self._get_results_array(grid_list)
        
            
        else:
            comm_grid = Grid(origin=self.common_grid_origin,
                             far_corner=self.common_grid_far_corner, 
                             spacing=0.5,
                             _grid = None)
            self.common_grid_nsteps = comm_grid.nsteps
                                                       
            results_array = np.zeros(self.common_grid_nsteps, dtype=tuple)
            assert len(origin_diff) == len(grid_list)
           
            for i in range(len(origin_diff)):
                g = grid_list[i]
                diff = origin_diff[i]
                d_x, d_y, d_z = diff
                nx, ny, nz = g.nsteps
                for x in range(nx):
                    for y in range(ny):
                        for z in range(nz):
                            if g.value(x, y, z) != 0:
                                if isinstance(results_array[x+ int(d_x)][y + int(d_y)][z + int(d_z)], tuple):
                                    results_array[x+ int(d_x)][y + int(d_y)][z + int(d_z)] += (g.value(x, y, z),)
                                else:
                                    results_array[x+ int(d_x)][y + int(d_y)][z + int(d_z)] = (g.value(x, y, z),)

            self.results_array = results_array
            self.nonzeros = self.results_array.nonzero()
                
                
    #### Functions for analysing ensemble data #####

    def get_gridpoint_means(self):
        """
        For each tuple in the GridEnsemble, calculates the difference in score between each point in the tuple and the mean of the tuple. 
        :return: list
        """
        means = [np.mean(val) for val in self.results_array[self.nonzeros]]
        return means

    def get_gridpoint_max(self):
        """
        For each tuple in the GridEnsemble, calculates the max of the tuple
        :return: list
        """
        maxes = [np.max(val) for val in self.results_array[self.nonzeros]]
        return maxes

    def get_gridpoint_ranges(self):
        """
        For each tuple in the GridEnsemble, returns the difference between max and mean
        :return: list
        """
        ranges = [np.max(val) - np.min(val) for val in self.results_array[self.nonzeros]]
        return ranges

    #### Functions for plotting histograms of analysed ensemble data ####

    def get_gridpoint_histograms(self):
        """
        Makes and saves histograms for each point in the results array.
        Caution - may output thousands of histograms for large maps.
        """

        ind_array = np.indices(self.results_array.shape)

        def results_array_histograms(x, y, z):
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Num_zeros: ', num_zeros)
                hist_arr = np.array(self.results_array[x][y][z])
                # hist, bin_edges = np.histogram(hist_arr, bins=20)
                colour_dict = {"acceptor": "r", "donor": "b", "apolar": "y"}
                hist_name = self.prot_name + '_' + self.probe + '_{}_{}_{}'.format(x, y, z)

                plt.figure(1)
                plt.hist(hist_arr, bins=20, color=colour_dict[self.probe])
                plt.figtext(0.6, 0.8, ('Number of zero values:' + str(num_zeros)))
                plt.title('Score distribution at point x:{}, y:{}, z:{}'.format(x, y, z))
                plt.xlabel('Fragment hotspot score')
                plt.ylabel('Frequency')
                plt.savefig(join(self.out_dir, hist_name))
                plt.close()

        print('Generating Histograms')
        vresults_array_histograms = np.vectorize(results_array_histograms)
        vresults_array_histograms(ind_array[0], ind_array[1], ind_array[2])


    def plot_gridpoint_spread(self):
        '''
        For each point in the 3D grid, plots the difference in score between each point in the tuple and the mean of the tuple. 

        '''
        means = self.get_gridpoint_means()
        mean_arr = np.array(means)
        (mu, sigma) = norm.fit(mean_arr)
        n, bins, patches = plt.hist(means, bins=40, normed=1)
        # print(bins)
        # y_fit = np.random.normal(mu, sigma, np.shape(mean_arr))
        y = norm.pdf(bins, mu, sigma)
        a = plt.plot(bins, y, 'r--', linewidth=2)
        bins = 0.5 * (bins[1:] + bins[:-1])
        y_fit = norm.pdf(bins, mu, sigma)
        ss_res = np.sum((n - y_fit) ** 2)
        ss_tot = np.sum((n - np.mean(n)) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        plt.title('Mu: {}, Sigma: {}, R^2 : {} '.format(round(mu, 2), round(sigma, 2), round(r2, 4)))
        hist_name = self.prot_name + '_{}_gridpoint_spread_fit'.format(self.probe)
        plt.savefig(join(self.out_dir, hist_name))
        # plt.show()
        plt.close()

    def plot_gridpoint_ranges(self):
        '''
        Plots the range of the tuple values for each point in the 3D grid
        :return: ranges = list of the tuple ranges at each point
        '''
        ranges = self.get_gridpoint_ranges()
        plt.hist(ranges, bins=40, normed=0)
        plt.title('Score ranges for {} {}'.format(self.prot_name, self.probe))
        hist_name = self.prot_name + '_{}_score_ranges'.format(self.probe)
        plt.savefig(join(self.out_dir, hist_name))
        # plt.show()
        plt.close()

    #### Functions that output analysed ensemble data as Grids ####

    def _make_grid(self, values):
        """
        Makes a grid to store output of ranges, max, means, etc 
        Duplicating Grid.from_array
        :param values: 
        :return: 
        """
        grid = Grid(origin=self.common_grid_origin,
                    far_corner=self.common_grid_far_corner,
                    spacing= 0.5,
                    default=0.0,
                    _grid=None)

        as_triads = zip(*self.nonzeros)
        for (i, j, k), v in zip(as_triads, values):
            grid._grid.set_value(int(i), int(j), int(k), v)

        return grid

    def output_grid(self, mode="max", save=True):
        """
        Ouputs ensemble information as ccdc.utilites.Grid 
        :param mode: the operation to be done on each tuple in the grid. Can be "max", "mean", "ranges", or "frequency" (Length of tuple)
        :param save: bool, whether output should be written to disk
        :return: a class ccdc.utilities.Grid object
        """
        # Get the values
        if mode == "max":
            vals = self.get_gridpoint_max()
        elif mode == "mean":
            vals = self.get_gridpoint_means()
        elif mode =="ranges":
            vals = self.get_gridpoint_ranges()
        elif mode == "frequency":
            vals = [np.max(val) - np.min(val) for val in self.results_array[self.nonzeros]]
        else:
            print("Unrecognised mode: {}".format(mode))
            return

        # Fill the grid
        out_grid = self._make_grid(vals)

        if save:
            out_grid.write(join(self.out_dir, '{}_{}_{}.ccp4'.format(self.prot_name,mode, self.probe)))

        return out_grid

    ###### Saving and loading ensembles #########

    @staticmethod
    def load_GridEnsemble(filename):
        """
        Loads a pickled MegaGrid
        :param filename: str, full path to pickled grid
        :return: MegaGrid object
        """
        pickle_file = open(filename, 'rb')
        newGridEnsemble = pickle.load(pickle_file)
        return newGridEnsemble

    def pickle_GridEnsemble(self):
        """
        Saves MegaGrids as pickles.
        :return: 
        """
        pickle.dump(self, open(join(self.out_dir, '{}_{}_GridEnsemble.p'.format(self.prot_name, self.probe)), 'wb'))

    ###### Functions to run full ensemble calculation and output grids (to integrate into main hotspots code ######

    def from_hotspot_maps(self, path_list, out_dir, prot_name, probe_name, mode="max"):
        """
        Creates a GridEnsemble from paths to Hotspot maps for a certain probe
        :param stem: path to directory where the grids are
        :param out_dir: path to where grids and histograms are saved
        :param prot_name: str
        :param probe_name: 'donor', 'acceptor', or 'apolar'
        :return: 
        """
        print('In from_hotspot_maps')
        self.path_list = path_list
        self.out_dir = out_dir
        self.prot_name = prot_name
        self.probe = probe_name

        #common_grids = self._common_grids_from_paths()
        grid_list = [Grid.from_file(p) for p in self.paths_list]
        self.get_alternative_results_array(grid_list)

        return self.output_grid(mode, save=True)

    def from_grid_list(self, grid_list, out_dir, prot_name, probe_name, mode="max"):
        """
        Creates a GridEnsemble from Hotspot maps for a certain probe
        :param grid_list: 
        :param out_dir: 
        :param prot_name: 
        :param probe_name: 
        :return: 
        """
        self.out_dir = out_dir
        self.prot_name = prot_name
        self.probe = probe_name
        print("Making ensemble {} {}".format(self.prot_name, self.probe))

        #common_grids = self._common_grids_from_grid_list(grid_list)
        #self._get_results_array(common_grids)
        self.get_alternative_results_array(grid_list)
        print(self.common_grid_origin, self.common_grid_far_corner)

        return self.output_grid(mode, save=False)

