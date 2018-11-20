'''
The main class of the :mod:`hotspots.grid_extension.Grid`.

This is an internal extension of :class:`ccdc.grid.Grid` that adds in potential new features
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
from hotspots.hotspot_utilities import Helper
from scipy import ndimage

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
    def from_molecule(mol, padding=1, scaling=0.5):
        """
        generate a molecule mask where gp within the vdw radius of the molecule heavy atoms are set to 1.0
        :param mol: `ccdc.molecule.Molecule`
        :param padding: int
        :param scaling: float
        :return: `hotspots.grid_extension.Grid`
        """
        x = set()
        y = set()
        z = set()
        for a in mol.atoms:
            x.add(a.coordinates.x)
            y.add(a.coordinates.y)
            z.add(a.coordinates.z)
        origin = Coordinates(x=round(min(x) - padding),
                             y=round(min(y) - padding),
                             z=round(min(z) - padding))

        far_corner = Coordinates(x=round(max(x) + padding),
                                 y=round(max(y) + padding),
                                 z=round(max(z) + padding))

        g = Grid(origin=origin, far_corner=far_corner, spacing=0.5, default=0, _grid=None)
        for a in mol.heavy_atoms:
            g.set_sphere(point=a.coordinates,
                         radius=a.vdw_radius * scaling,
                         value=1,
                         scaling='None')
        return g > 0.1

    @staticmethod
    def grow(inner, template, percentile=60):
        """
        experimental
        Dilates grid to the points in the top percentile of the template
        :param template:
        :return:
        """
        expand = inner.max_value_of_neighbours() > 0.2
        outer = expand.__sub__(inner) * template
        threshold = np.percentile(a=outer.grid_values(threshold=1), q=int(percentile))
        return inner.__add__(outer > threshold)

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
