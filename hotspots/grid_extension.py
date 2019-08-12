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
from skimage.morphology import ball
from os.path import join, basename
from scipy.stats import norm
import matplotlib.pyplot as plt
import pickle
from functools import reduce

Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


class Grid(utilities.Grid):
    """
    A class to extend a `ccdc.utilities.Grid` this provides grid methods required in the Fragment Hotspot Maps algorithm
    """
    def coordinates(self, threshold=1):
        nx, ny, nz = self.nsteps
        return [self.indices_to_point(i, j, k)
                for i in range(nx)
                for j in range(ny)
                for k in range(nz)
                if self.value(i, j, k) >= threshold]


    def grid_value_by_coordinates(self, threshold=1):
        """
        returns a dictionary of grid point values by coordinates
        :param int threshold: the island threshold of the grid, only points over this value will be returned
        :return: dict, grid point values by coordinates
        """
        dic = {}
        nx, ny, nz = self.nsteps
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if self.value(i, j, k) > threshold:
                        try:
                            dic[self.value(i, j, k)].append(self.indices_to_point(i, j, k))
                        except KeyError:
                            dic.update({self.value(i, j, k): [self.indices_to_point(i, j, k)]})

        return dic

    @staticmethod
    def _tolerance_range(value, tolerance, min, max):
        """
        private method

        for a given value and tolerance, the method checks that the tolerance range for that value is within the grid
        boundaries
        :param int value: an "indice" value either (i or j or k)
        :param int min: the minimum grid boundary (always = 0)
        :param int max: the maximum grid boundary (always = n step)
        :return: range,
        """
        low = value - tolerance
        if low < 0:
            low = 0

        high = value + tolerance
        if high > max:
            high = max
        return range(low, high)

    def get_near_scores(self, coordinate, tolerance=3):
        """
        for a given grid point, return a list of values within a search tolerance
        :param tup coordinate: coordinate of a point within the grid
        :param int tolerance: search distance, in grid steps
        :return:
        """
        i, j, k = self.point_to_indices(coordinate)
        ri = self._tolerance_range(i, tolerance, 0, self.nsteps[0])
        rj = self._tolerance_range(j, tolerance, 0, self.nsteps[1])
        rk = self._tolerance_range(k, tolerance, 0, self.nsteps[2])
        return [self.value(a, b, c) for a in ri for b in rj for c in rk if self.value(a,b,c) > 0]

    def grid_values(self, threshold=0):
        """
        generates a numpy array with all values over a given threshold (default = 0)
        :param int threshold: values over this value
        :return:
        """
        array = self.get_array()
        masked_array = np.ma.masked_less_equal(array, threshold)
        return masked_array.compressed()

    def grid_score(self, threshold=0, percentile=75):
        """
        for a given grid, the xth percentile of values above a given threshold is returned
        :param int threshold: values over this value
        :param int percentile: value at this percentile

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
        return x,y,z coordinate for a given grid index
        :param int i: indice (x-axis)
        :param int j: indice (y-axis)
        :param int k: indice (z-axis)
        :return: tup, float(x), float(y), float(z)
        """
        ox, oy, oz, = self.bounding_box[0]
        gs = self.spacing
        return ox + float(i) * gs, oy + float(j) * gs, oz + gs * float(k)

    def point_to_indices(self, p):
        """
        return the nearest grid index for a given point
        :param tup p: (float(x), float(y), float(z)), a coordinate on a 3D grid
        :return: int(x), int(y), int(z)
        """

        gs = self.spacing
        rx, ry, rz = [round(i / gs) for i in p]
        ox, oy, oz = [round(i / gs) for i in self.bounding_box[0]]
        return int(rx - ox), int(ry - oy), int(rz - oz)

    def gaussian(self, sigma=0.2):
        """
        gaussian smoothing function, method of reducing noise in output
        :param float sigma: degree of smoothing
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
        :param tup point: (float(x), float(y), float(z))
        :param int threshold: values above this value
        :param float tolerance: radius of search
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
        :return: `numpy.array`, array with shape (nsteps) and each element corresponding to value at that indice
        """
        nx, ny, nz = self.nsteps
        array = np.zeros((nx, ny, nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    array[i, j, k] += self.value(i, j, k)
        return array

    def dilate_by_atom(self):

        g_array = self.get_array()
        selem = ball(radius=2)
        print(selem)

        dilated = ndimage.grey_dilation(g_array, structure=selem)

        return self.array_to_grid(dilated,self)




    def restricted_volume(self, volume=75):
        """
        returns a grid with of a defined volume
        :param float volume: desired volume in Angstroms ^ 3
        :return: `hotspots.grid_extension.Grid`
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
        returns centre of the grid's bounding box
        :param self:
        :return:
        """
        return ((self.bounding_box[0][0] + self.bounding_box[1][0])/2,
                (self.bounding_box[0][1] + self.bounding_box[1][1])/2,
                (self.bounding_box[0][2] + self.bounding_box[1][2])/2
                )

    def deduplicate(self, major, threshold=12, tolerance=2):
        """
        method to deduplicate two grids, used for charged-polar deduplication
        :param `ccdc.utilities.Grid` major: overriding grid
        :param int threshold: values above this value
        :param int tolerance: search radius for determining feature overlap
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
        make a new empty grid, with the same dimensions as the old
        :return: `hotspots.grid_extension.Grid`
        """
        g = self.copy()
        g *= 0
        return g

    def common_boundaries(self, grid):
        """
        expands supplied grid to the size of self (supplied grid should be smaller than self)
        :param grid:
        :return:
        """
        reference = Grid.super_grid(0, self, grid)
        blank = reference.copy_and_clear()
        return Grid.super_grid(0, blank, grid)

    def multi_max_mask(self, grids):
        """
        for a self grid and collection of grids supplied, for each grid point, if the maximum grid value across all the
        grids belongs to the self grid, the value is assigned to a fresh grid.
        :param list grids: grids to be compared to self
        :return: `ccdc.utilities.Grid`
        """
        max_grids = [self > g for g in grids]
        blank = -self.copy_and_clear()
        return reduce(operator.__and__, max_grids, blank)

    def get_best_island(self, threshold, mode="count", peak=None):
        """
        returns the best grid island. Mode: "count" or "score"
        :param int threshold: island threshold
        :param str mode:
                       -"count" : returns island with most grid points above threshold
                       -"score" : returns island with the largest sum of all grid points over threshold
        :param tup peak: (float(x), float(y), float(z)) coordinates of peak in grid
        :return: `ccdc.utilities.Grid`, grid containing the best island
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
                            island_by_rank.update({rank: g})
                        else:
                            continue
                    else:
                        g = (island > threshold) * island
                        rank = g.count_grid()
                        island_by_rank.update({rank: g})

            # elif mode == "score":
            #     for island in islands:
            #         if peak:
            #             if island.contains_point(peak):
            #                 nx, ny, nz = island.nsteps
            #                 island_points = [island.value(i, j, k)
            #                                  for i in range(nx) for j in range(ny) for k in range(nz)
            #                                  if island.value(i, j, k) >= threshold]
            #                 rank = sum(island_points)
            #                 island_by_rank.update({rank: island})
            #             else:
            #                 continue
            #         else:
            #             nx, ny, nz = island.nsteps
            #             island_points = [island.value(i, j, k)
            #                              for i in range(nx) for j in range(ny) for k in range(nz)
            #                              if island.value(i, j, k) >= threshold]
            #             rank = sum(island_points)
            #             island_by_rank.update({rank: island})

            else:
                raise IOError("mode not supported")

            if len(island_by_rank) == 0:
                return None
            else:
                rank = sorted(island_by_rank.keys(), reverse=True)[0]
                print("threshold:", threshold, "count:", sorted(island_by_rank.keys(), reverse=True))
                return island_by_rank[rank]

    def minimal(self):
        """
        reduces grid size to the minimal dimensions
        :return: `ccdc.utilities.Grid`
        """
        try:
            return Grid.super_grid(1, *self.islands(threshold=1))
        except RuntimeError:
            return self

    def limit_island_size(self, npoints, threshold=10):
        """
        for a given grid, the number of points contained within the islands (threshold = x) is limited to npoints
        :param int npoints: maximum number of points in each island feature
        :param float threshold: values above this value
        :return:
        """
        g = (self > 10) * self
        all_islands = []
        for island in g.islands(threshold):
            if island.count_grid() > npoints:
                all_islands.append(island.top_points(npoints=npoints))
            else:
                all_islands.append(island)
        return Grid.super_grid(0, *all_islands)

    def top_points(self, npoints):
        """
        for a given grid, the top scoring n points are returned
        :param int npoints: number of points to be returned
        :return: `ccdc.ulilities.Grid`
        """
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
        return (self > thres) * self

    def step_out_mask(self, nsteps=2):
        """
        add n step in all directions to Grid boundary, returns blank grid
        :param nsteps: the number of steps the grid is expanded
        :return: `hotspots.grid_extension.Grid`
        """
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
        """
        private method
        returns the lower limit of the top n number of points in a grid
        :param dict sorted_points: {grid_value, (float(x), float(y), float(z))
        :param int npoints: number of grid points
        :return:
        """
        count = []
        for value, pts in sorted_points:
            count.extend(pts)
            if len(count) >= npoints:
                return value
            else:
                continue

    @staticmethod
    def array_to_grid(array, blank):
        """"""
        grid = blank.copy_and_clear()
        indices = np.nonzero(array)
        values = array[indices]
        as_triads = zip(*indices)

        for (i, j, k), v in zip(as_triads, values):

            grid._grid.set_value(int(i), int(j), int(k), float(v))

        return grid

    @staticmethod
    def from_array(fname):
        """
        creates a grid from array
        :param fname: path to pickled numpy array
        :return: `hotspots.grid_extension.Grid`
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
        :param list grid_list: list of `ccdc.utilities.Grid` instances
        :param padding: number of steps to add to the grid boundary
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
            return mask_dic, reduce(operator.add, mask_dic.values(), blank).minimal()
        else:
            return reduce(operator.add, mask_dic.values(), blank).minimal()

    def _mutually_inclusive(self, other):
        """

        :param other:
        :return:
        """
        g, h = Grid.common_grid(grid_list=[self, other], padding=1)
        return g & h

    def atomic_overlap(self, atom, return_grid=True):

        g = Grid.initalise_grid(coords=[atom.coordinates])

        h = g.copy_and_clear()
        h.set_sphere(point=atom.coordinates, radius=atom.vdw_radius, value=1, scaling='None')
        overlap = self._mutually_inclusive(other=h)
        perc_overlap = (overlap.count_grid() / (h > 0).count_grid()) * 100

        if return_grid is True:
            return perc_overlap, overlap

        else:
            return perc_overlap

    def matched_atoms(self, atoms, threshold=30):
        """
        for a given atom, the percentage overlap with the grid is calculated. If the overlap
        is over a threshold the atom identifier is returned in a list

        :param list atoms: list of `ccdc.molecule.Atoms`
        :param int threshold: percentage overlap threshold
        :return list: list of str
        """
        passed_atoms = {}

        for a in atoms:
            perc_overlap, overlap = self.atomic_overlap(atom=a, return_grid=True)

            if perc_overlap > threshold:
                common_a, common_b = Grid.common_grid([self, overlap])
                passed_atoms[a.label] = (common_a * common_b).extrema[1]

        return passed_atoms

    def percentage_overlap(self, other):
        """
        find the percentage overlap of this grid with other.
        :param other: `hotspots.grid_extension.Grid`
        :return:`hotspots.grid_extension.Grid`
        """
        overlap = self._mutually_inclusive(other=other).count_grid()
        vol = (self > 0).count_grid()
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
        coords = [a.coordinates for a in mol.atoms]
        g = Grid.initalise_grid(coords=coords)
        for a in mol.heavy_atoms:
            g.set_sphere(point=a.coordinates,
                         radius=a.vdw_radius * scaling,
                         value=1,
                         scaling='None')
        return g > 0.1

    @staticmethod
    def initalise_grid(coords, padding=1, spacing=0.5):
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

        return Grid(origin=origin, far_corner=far_corner, spacing=spacing, default=0, _grid=None)

    @staticmethod
    def grow(inner, template, percentile= 80):
        """
        experimental
        Dilates grid to the points in the top percentile of the template
        :param template:
        :return:
        """
        expand = inner.max_value_of_neighbours() > 0.1   # remove very small values
        outer = expand.__mul__(-inner) * template
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
    """
    Experimental feature
    Class that handles a numpy array of tuples from Hotspot maps of a given probe type across multiple structures
    """

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
        self.tup_max_length = len(grid_list)
        self.common_grid_origin = common_grids[0].bounding_box[0]
        self.common_grid_far_corner = common_grids[0].bounding_box[1]
        self.common_grid_nsteps = common_grids[0].nsteps


        return common_grids


    def get_4D_results_array(self, grid_list):
        """
         Reads in grids, converts them to 3d numpy arrays, and stacks them into 4d numpy array, which holds the information
         for the ensemble.
         :return:
         """
        # Initialise the array
        common_grids = self._common_grids_from_grid_list(grid_list)
        common_arrays = [cg.get_array() for cg in common_grids]
        ensemble_array = np.stack(common_arrays, axis=3)
        self.results_array = ensemble_array

        max_array = np.max(self.results_array, axis=-1)
        self.nonzeros = max_array.nonzero()


    #### Functions for analysing ensemble data #####

    def get_gridpoint_means(self):
        """
        For each tuple in the GridEnsemble, calculates the means of the tuple
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

    def get_gridpoint_means_spread(self):
        """
        For tuple in the GridEnsemble, calculates the difference in score between each point in the tuple and the mean of the tuple.
        :return: Python list
        """
        means_spread = []
        values = self.results_array[self.nonzeros]

        for v in values:
            num_zeros = self.tup_max_length - len(v)
            if num_zeros != 0:
                print('Number of zeros', num_zeros, 'out of ', self.tup_max_length)
            hist_arr = np.array(v)
            means_spread.extend(list(hist_arr - np.mean(hist_arr)))

        return means_spread

    # Functions for plotting histograms of analysed ensemble data ####
    def plot_gridpoint_spread(self):
        '''
        For each point in the 3D grid, plots the difference in score between each point in the tuple and the mean of the tuple.
        '''
        # Get the data:
        means_spread = self.get_gridpoint_means_spread()
        mean_arr = np.array(means_spread)
        (mu, sigma) = norm.fit(mean_arr)
        n, bins, patches = plt.hist(means_spread, bins=40, normed=1)
        # print(bins)
        # y_fit = np.random.normal(mu, sigma, np.shape(mean_arr))
        y = norm.pdf(bins, mu, sigma)
        plt.plot(bins, y, 'r--', linewidth=2)
        fit_bins = 0.5 * (bins[1:] + bins[:-1])
        y_fit = norm.pdf(fit_bins, mu, sigma)
        ss_res = np.sum((n - y_fit) ** 2)
        ss_tot = np.sum((n - np.mean(n)) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        plt.title('Mu: {}, Sigma: {}, R^2 : {} '.format(round(mu, 2), round(sigma, 2), round(r2, 2)))
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
                    spacing=0.5,
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
        elif mode == "ranges":
            vals = self.get_gridpoint_ranges()
        elif mode == "frequency":
            vals = [np.max(val) - np.min(val) for val in self.results_array[self.nonzeros]]
        else:
            print("Unrecognised mode: {}".format(mode))
            return

        # Fill the grid
        out_grid = self._make_grid(vals)

        if save:
            out_grid.write(join(self.out_dir, '{}_{}_{}.ccp4'.format(self.prot_name, mode, self.probe)))

        return out_grid

    ###### Saving and loading ensembles #########

    @staticmethod
    def load_GridEnsemble(filename):
        """
        Loads a pickled _GridEnsemble
        :param filename: str, full path to pickled grid
        :return: MegaGrid object
        """
        pickle_file = open(filename, 'rb')
        newGridEnsemble = pickle.load(pickle_file)
        return newGridEnsemble

    def pickle_GridEnsemble(self):
        """
        Saves _GridEnsembles as pickles.
        :return:
        """
        pickle.dump(self, open(join(self.out_dir, '{}_{}_GridEnsemble.p'.format(self.prot_name, self.probe)), 'wb'))

    ###### Functions to run full ensemble calculation and output grids (to integrate into main hotspots code ######

    def from_hotspot_maps(self, path_list, out_dir, prot_name, probe_name, mode="max"):
        """
        Creates a GridEnsemble from paths to Hotspot maps for a certain probe
        :param path_list: list of paths for the hotspot maps
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

        grid_list = [Grid.from_file(p) for p in self.path_list]
        # self.get_alternative_results_array(grid_list)
        self.get_4D_results_array(grid_list)

        return self.output_grid(mode, save=False)

    def from_grid_list(self, grid_list, out_dir, prot_name, probe_name, mode="max"):
        """
        Creates a GridEnsemble from Hotspot maps for a certain probe
        :param grid_list:
        :param out_dir:
        :param prot_name:
        :param probe_name:
        :return:
        """
        print(mode)
        self.out_dir = out_dir
        self.prot_name = prot_name
        self.probe = probe_name
        print("Making ensemble {} {}".format(self.prot_name, self.probe))

        # common_grids = self._common_grids_from_grid_list(grid_list)
        self.get_4D_results_array(grid_list)
        print(self.common_grid_origin, self.common_grid_far_corner)

        return self.output_grid(mode, save=False)

