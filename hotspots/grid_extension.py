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
from scipy.spatial import distance
from skimage import feature
from skimage.morphology import ball
from os.path import join, basename
from scipy.stats import norm
import matplotlib.pyplot as plt
import pickle
from functools import reduce

from skimage.morphology import remove_small_objects

try:
    from hdbscan import HDBSCAN
except ImportError:
    print('HDBSCAN module not installed. _GridEnsemble clustering not available.')

Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


class Grid(utilities.Grid):
    """
    A class to extend a `ccdc.utilities.Grid` this provides grid methods required in the Fragment Hotspot Maps algorithm
    """

    # creating new grids
    ##################################################################################################################
    @staticmethod
    def from_molecule(mol, scaling=1, value=1, scaling_type='None', spacing=0.5):
        """
        generate a molecule mask where gp within the vdw radius of the molecule heavy atoms are set to 1.0
        :param mol: a molecule
        :type mol: :class:`ccdc.molecule.Molecule`
        :param scaling: scale radius by this value
        :type scaling: float
        :param value: value to set grid point to.
        :type value: float
        :param scaling_type: scalling of grid point values from centre to edge. Either 'none' or 'linear'
        :type scaling_type: str
        :return: a grid
        :rtype: :class:`hotspots.grid_extension.Grid`
        """
        coords = [a.coordinates for a in mol.atoms]
        g = Grid.initalise_grid(coords=coords, padding=2, spacing=spacing)
        for a in mol.heavy_atoms:
            g.set_sphere(point=a.coordinates,
                         radius=a.vdw_radius * scaling,
                         value=value,
                         scaling=scaling_type)
        return g

    def from_coords(self, coords, scaling=1):
        """
        generate a molecule mask where gp within the vdw radius of the molecule heavy atoms are set to 1.0
        :param mol: `ccdc.molecule.Molecule`
        :param padding: int
        :param scaling: float
        :return: `hotspots.grid_extension.Grid`
        """

        g = self.copy()
        for c in coords:
            g.set_sphere(point=c,
                         radius=1 * scaling,
                         value=1,
                         scaling='None')
        return g

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
        try:
            for c in coords:
                x.add(c.x)
                y.add(c.y)
                z.add(c.z)
        except:
            for c in coords:
                x.add(c[0])
                y.add(c[1])
                z.add(c[2])

        origin = Coordinates(x=round(min(x) - padding),
                             y=round(min(y) - padding),
                             z=round(min(z) - padding))

        far_corner = Coordinates(x=round(max(x) + padding),
                                 y=round(max(y) + padding),
                                 z=round(max(z) + padding))

        return Grid(origin=origin, far_corner=far_corner, spacing=spacing, default=0, _grid=None)

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




    ##################################################################################################################
    # reformat grid data: methods which restructure grid data i.e. grid to numpy array
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

    def grid_values(self, threshold=0):
        """
        generates a numpy array with all values over a given threshold (default = 0)
        :param int threshold: values over this value
        :return:
        """
        array = self.get_array()
        masked_array = np.ma.masked_less_equal(array, threshold)
        return masked_array.compressed()

    def get_flat_array(self):
        """

        :return:
        """
        return np.array(self.to_vector())

    def get_array(self):
        """
        convert grid object to np.array
        :return: `numpy.array`, array with shape (nsteps) and each element corresponding to value at that indice
        """
        nx, ny, nz = self.nsteps
        flat = np.array(self.to_vector())
        array = flat.reshape((nx,ny,nz))
        return array

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
    def from_array(fname, origin, spacing):
        """
        creates a grid from array
        :param fname: path to pickled numpy array
        :return: `hotspots.grid_extension.Grid`
        """
        array = np.load(fname)
        shape = array.shape
        far_corner = [origin[0]+(spacing*shape[0]), origin[1]+(spacing*shape[1]), origin[2]+(spacing*shape[2])]

        grid = Grid(origin=list(origin),
                    far_corner=list(far_corner),
                    spacing=spacing,
                    default=0.0,
                    _grid=None)

        indices = np.nonzero(array)
        values = array[indices]
        as_triads = zip(*indices)

        for (i, j, k), v in zip(as_triads, values):
            grid._grid.set_value(int(i), int(j), int(k), float(v))

        return grid

    ##################################################################################################################
    # manipulate grid dimension: everything concerned with altering dimensions of a grid
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

    def minimal(self):
        """
        reduces grid size to the minimal dimensions
        :return: `ccdc.utilities.Grid`
        """
        try:
            return Grid.super_grid(1, *self.islands(threshold=1))
        except RuntimeError:
            return self

    @staticmethod
    def shrink(small, big, reverse_padding=1):
        """
        shrink a big grid to the dimension of a small grid

        :param big: the grid to be shrunk
        :type big: `hotspots.grid_extension.Grid`
        :param reverse_padding: amount of erosion within the small grid boundaries (ensures fit preventing a seg fault)
        :type reverse_padding: int

        :return: shrunk grid
        :rtype: `hotspots.grid_extension.Grid`
        """

        origin, far_left = small.bounding_box
        o = big.point_to_indices(origin)
        o = [i + reverse_padding for i in o]

        f = big.point_to_indices(far_left)
        f = [i - reverse_padding for i in f]

        h = big.sub_grid(o + f)
        # reverse padding ensure h is smaller than 'small'. Finally expand h to the dimensions of small.
        return small.common_boundaries(h)




    ##################################################################################################################
    # grid data calls: methods which access a part of grid data
    def centroid(self):
        """
        returns centre of the grid's bounding box
        :param self:
        :return:
        """
        return ((self.bounding_box[0][0] + self.bounding_box[1][0]) / 2,
                (self.bounding_box[0][1] + self.bounding_box[1][1]) / 2,
                (self.bounding_box[0][2] + self.bounding_box[1][2]) / 2
                )

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
        return [self.value(a, b, c) for a in ri for b in rj for c in rk if self.value(a, b, c) > 0]

    @staticmethod
    def neighbourhood(i, j, k, high, catchment=1):
        """
        find the neighbourhood of a given indice. Neighbourhood is defined by all points within 1 step of the
        specified indice. This includes the cubic diagonals.

        :param i: i indice
        :param j: j indice
        :param k: k indice
        :param catchment: number of steps from the centre

        :type i: int
        :type j: int
        :type k: int
        :type catchment: int

        :return: indices of the neighbourhood
        :rtype: list
        """
        low = (0, 0, 0)

        i_values = [a for a in range(i-catchment, i+catchment+1) if low[0] <= a < high[0]]
        j_values = [b for b in range(j-catchment, j+catchment+1) if low[1] <= b < high[1]]
        k_values = [c for c in range(k-catchment, k+catchment+1) if low[2] <= c < high[2]]

        return [[a, b, c] for a in i_values for b in j_values for c in k_values
                if Helper.get_distance([a, b, c], [i, j, k]) == 1]

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

    ##################################################################################################################
    # grid data manipulation
    @staticmethod
    def grow(inner, template, percentile=80):
        """
        experimental
        Dilates grid to the points in the top percentile of the template
        :param template:
        :return:
        """
        expand = inner.max_value_of_neighbours() > 0.1  # remove very small values
        outer = expand.__mul__(-inner) * template
        threshold = np.percentile(a=outer.grid_values(threshold=1), q=int(percentile))

        return inner.__add__(outer > threshold)

    def get_peaks(self, min_distance=6, cutoff=2):
        """
        -     Local maxima with at least a seperation of ((2 * min_distance) + 1)
              are returned
        -     If there are multiple local maxima with identical pixel intensities inside
              the region defined with `min_distance`, the coordinates of all such pixels
              are returned.
        -     Therefore, local maxima with multiple pixels can be grouped by distance
        :return:
        """
        class Peak:
            def __init__(self, score, indices):
                self.score = score
                self.indices = [indices]

            def centroid(self):
                x = set()
                y = set()
                z = set()

                for i in self.indices:
                    x.add(i[0])
                    y.add(i[1])
                    z.add(i[2])
                return [sum(x) / len(x), sum(y) / len(y), sum(z) / len(z)]

        peaks = feature.peak_local_max(self.get_array(), min_distance=min_distance, threshold_abs=cutoff)

        grouped_peaks = []
        threshold = (2 * min_distance) + 1

        for i, peak in enumerate(peaks):
            x, y, z = peak

            if i == 0:
                grouped_peaks.append(Peak(score=self.value(int(x), int(y), int(z)), indices=peak))

            else:

                min_d = [x < threshold for x in [np.amin(distance.cdist(np.array([peak]),
                                                                        np.array(g.indices)))
                                                 for g in grouped_peaks]
                         ]

                if any(min_d):
                    loci = (np.array(min_d) * 1).nonzero()
                    if len(loci) == 1:
                        x = loci[0][0]
                    else:
                        raise NotImplemented
                    grouped_peaks[x].indices.append(peak)

                else:
                    grouped_peaks.append(Peak(score=self.value(int(x), int(y), int(z)), indices=peak))

        average_peaks = []
        for p in grouped_peaks:
            i, j, k = p.centroid()
            coords = self.indices_to_point(i, j, k)
            average_peaks.append(coords)

        return average_peaks

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

        spacing = grd_dict["apolar"]._spacing

        blank = Grid(origin=origin, far_corner=far_corner, spacing=spacing, default=0.1, _grid=None)

        if mask:
            return mask_dic, reduce(operator.add, mask_dic.values(), blank)
        else:
            return reduce(operator.add, mask_dic.values(), blank)

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

    def dilate_by_atom(self, radius=1):

        g_array = self.get_array()
        selem = ball(radius=radius)
        print(selem)

        dilated = ndimage.grey_dilation(g_array, structure=selem)

        return self.array_to_grid(dilated, self)

    def get_best_island(self, threshold):
        """
        For a given threshold, the island which contains the most grid points will be returned

        :param threshold: island threshold
        :type threshold: int

        :return: the island with the most grid points above the threshold
        :rtype: :class:`ccdc.utilities.Grid`
        """
        islands = self.islands(threshold)

        if len(islands) == 0:
            return None

        else:
            island_by_rank = {}
            for island in islands:
                g = (island > threshold) * island
                rank = g.count_grid()
                island_by_rank.update({rank: g})

            if len(island_by_rank) == 0:
                return None

            else:
                rank = sorted(island_by_rank.keys(), reverse=True)[0]
                print("threshold:", threshold, "count:", sorted(island_by_rank.keys(), reverse=True))
                return island_by_rank[rank]

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

    def _mutually_inclusive(self, other):
        """

        :param other:
        :return:
        """
        g, h = Grid.common_grid(grid_list=[self, other], padding=1)
        return g & h

    def remove_small_objects(self, min_size = 400):

        g_array = self.get_array()
        bool_array = g_array.astype(bool)
        no_small_obj = remove_small_objects(bool_array, min_size=min_size, connectivity=2)

        out_array = g_array*no_small_obj

        return self.array_to_grid(out_array,self.copy_and_clear())

    ##################################################################################################################
    def score_atom(self, atom):
        selem = ball(radius=1)
        i,j,k = self.point_to_indices(atom.coordinates)
        score = 0
        len_i, len_j, len_k = np.shape(selem)

        for di in range(0,len_i):
            for dj in range(0,len_j):
                for dk in range(0,len_k):
                    # print(di, dj, dk)
                    if selem[di,dj,dk] ==1:

                        x = i+(di-2)
                        y = j+(dj-2)
                        z = k+(dj-2)
                        # print(x,y,z)
                        score += self.value(x,y,z)
        return score

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

    ##################################################################################################################
    # other grid
    def edge_detection(self, edge_definition=0):
        """
        A simplified method to detect surface edge. An edge is defined as a grid point has a value 1 but is adjacent
        to a grid point with value 0. Only points distance = 1 are considered adjacent (i.e. not diagonals)

        :param edge_definition: values above which are considered part of the body
        :type edge_definition: float

        :return: the bodies surface as a list of indices
        :rtype: list
        """
        edge = []
        # generate a mask.
        self = self > edge_definition
        a = self.get_array()
        nx, ny, nz = self.nsteps
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if a[i][j][k] > 0:
                        neighbourhood = self.neighbourhood(i, j, k, self.nsteps)
                        if min({a[n[0]][n[1]][n[2]] for n in neighbourhood}) == 0:
                            edge.append(self.indices_to_point(i, j, k))

        return edge

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

    # def restricted_volume(self, volume=75):
    #     """
    #     returns a grid with of a defined volume
    #     :param float volume: desired volume in Angstroms ^ 3
    #     :return: `hotspots.grid_extension.Grid`
    #     """
    #     grid = self.copy_and_clear()
    #     max_points = int(float(volume) / 0.125)
    #     nx, ny, nz = self.nsteps
    #     rank_dict = {}
    #
    #     for i in range(nx):
    #         for j in range(ny):
    #             for k in range(nz):
    #                 value = float(self.value(i, j, k))
    #                 if value in rank_dict:
    #                     rank_dict[value].append((i, j, k))
    #                 else:
    #                     rank_dict.update({value: [(i, j, k)]})
    #
    #     top_points = sorted((float(x) for x, y in rank_dict.iteritems()), reverse=True)
    #     indices = [pts for key in top_points for pts in rank_dict[key]]
    #
    #     for i in indices[:max_points]:
    #         grid.set_value(i[0], i[1], i[2], self.value(i[0], i[1], i[2]))
    #
    #     return Grid.super_grid(1, *grid.islands(threshold=1))

    # def limit_island_size(self, npoints, threshold=10):
    #     """
    #     for a given grid, the number of points contained within the islands (threshold = x) is limited to npoints
    #     :param int npoints: maximum number of points in each island feature
    #     :param float threshold: values above this value
    #     :return:
    #     """
    #     g = (self > 10) * self
    #     all_islands = []
    #     for island in g.islands(threshold):
    #         if island.count_grid() > npoints:
    #             all_islands.append(island.top_points(npoints=npoints))
    #         else:
    #             all_islands.append(island)
    #     return Grid.super_grid(0, *all_islands)


utilities.Grid = Grid


class _GridEnsemble(object):
    """
    Given a list of hotspot maps for the same probe type, compiles a 4-dimensional numpy array storing the information.
    """

    def __init__(self, dimensions=None, shape=None, ensemble_array=None, spacing=0.5):
        """

        :param grid_list: List of hotspot maps for a certain probe type
        :type list:
        """
        self.dimensions = dimensions
        self.shape = shape
        self.ensemble_array = ensemble_array
        self.spacing = spacing

    @staticmethod
    def array_from_grid(grid):
        """
        Converts a grid into a 3D nunmpy array

        :param grid:
        :return:
        """
        nx, ny, nz = grid.nsteps
        array = np.zeros((nx, ny, nz))
        grid_vec = grid.to_vector()

        array = np.array(grid_vec).reshape(grid.nsteps)
        return array

    def make_ensemble_array(self, grid_list):
        """
        Creates a 4D numpy array storing information for the ensemble

        :param grid_list:
        :return:
        """
        print("Making the common grids")
        common_grids = Grid.common_grid(grid_list)
        print("Stared making arrays")
        as_arrays = [self.array_from_grid(cg) for cg in common_grids]

        self.ensemble_array = np.stack(as_arrays, axis=-1)
        print("GridEnsemble complete")
        self.dimensions = np.array(common_grids[0].bounding_box)
        self.shape = common_grids[0].nsteps

    def as_grid(self, array):
        """
        Given an array, outputs a grid with the dimensions of the GridEnsemble

        :param array: 3D numpy array, usually containing processed ensemble data
        :return: a :class: 'ccdc.utilities.Grid' instance
        """
        # Initialise the Grid
        grid = Grid(origin=tuple(self.dimensions[0]),
                    far_corner=tuple(self.dimensions[1]),
                    spacing=self.spacing,
                    default=0.0,
                    _grid=None)
        # Get the nonzero indices and values of the array
        nonz = array.nonzero()
        values = array[nonz]
        # Get indices per value
        as_triads = zip(*nonz)

        # Fill in the grid
        for (i, j, k), v in zip(as_triads, values):
            grid._grid.set_value(int(i), int(j), int(k), v)

        return grid

    def make_summary_grid(self, mode='median'):
        """
        Returns a ccdc grid, containing information at each point according to the mode pamater

        :param mode: one of 'mean', 'max', or 'median' (default)
        :type str:

        :return: 'hotspots.grid_extension.Grid'
        """
        if mode == 'median':
            arr = np.median(self.ensemble_array, axis=-1)

        elif mode == 'mean':
            arr = np.mean(self.ensemble_array, axis=-1)

        elif mode == 'max':
            arr = np.max(self.ensemble_array, axis=-1)

        else:
            print('unrecognised mode of combining grids')
            return

        agrid = self.as_grid(arr)

        return agrid

    def make_nonzero_median_map(self):
        """
        Takes the median of only points >0 at each point in the ensemble
        :return: numpy array
        """
        nan_arr = self.ensemble_array.copy()
        nan_arr[nan_arr == 0] = np.nan
        med_map = np.nanmedian(nan_arr, axis=-1)
        med_map = np.nan_to_num(med_map)
        return med_map

    def get_frequency_map(self):
        """
        Calculates the percent of nonzero values at each point
        Equivalent of len(nonzero values)/len(all_values)*100 at every point
        :return: numpy array
        """
        return np.divide(np.count_nonzero(self.ensemble_array, axis=-1), float(self.ensemble_array.shape[-1])) * 100

    def get_difference_frequency_map(self, other, threshold):
        """

        :param other:
        :param threshold:
        :return: 3d numpy array
        """
        freq_map = self.get_frequency_map()
        med_diff_map = self.make_nonzero_median_map() - other.make_nonzero_median_map()
        thresh_med_map = med_diff_map * (freq_map > threshold)
        return thresh_med_map

    def get_median_frequency_map(self, threshold):
        """

        :param threshold:
        :return: 3d numpy array
        """
        freq_map = self.get_frequency_map()
        med_map = self.make_nonzero_median_map()
        thresh_med_map = med_map * (freq_map > threshold)

        return thresh_med_map

    @staticmethod
    def get_center_of_mass(array):
        """

        :param array:
        :return: coordinates as 1d numpy array
        """
        indices = array.nonzero()
        vals = array[indices]
        result_list = []
        for dimension in indices:
            result_list.append(np.average(dimension, weights=vals))
        return result_list

    @staticmethod
    def get_highest_percentile_scores(arr, percentile=90.0):
        """
        Returns a grid thresholded at the highest <percentile> of scores.
        :param percentile:
        :return:
        """
        vals = arr[arr.nonzero()]
        perc = np.percentile(vals, percentile)
        perc_arr = arr * (arr > perc)
        return perc_arr

    def get_contributing_maps(self, cluster_array):
        """
        Given an array with the same first 3 dimensions of the ensemble_array, with points labelled by cluster,
        returns a dictionary of which structures contribute to which cluster.
        :param cluster_array: 3D array, labelled by cluster (eg output of self.HDBSCAN_cluster())
        :return:
        """
        clust_dic = {}
        clusters = list(set(cluster_array[cluster_array.nonzero()]))

        for c in clusters:
            cluster_mask = (cluster_array == c)
            ensemble_cluster = self.ensemble_array[cluster_mask]  # should result in 2D array
            grid_indices = list(ensemble_cluster.nonzero()[1])
            clust_structs = list(set(grid_indices))
            clust_dic[c] = [(val, grid_indices.count(val)) for val in clust_structs]
        return clust_dic

    @staticmethod
    def HDBSCAN_cluster(d_array, **kwargs):
        """
        Performs density-based clustering on input 3d map.
        :param d_array: input numpy array (usually 3d map)
        :param min_members:
        :return: numpy array
        """
        clusterer = HDBSCAN(**kwargs)
        in_arr = np.array(d_array.nonzero()).T

        clusterer.fit(in_arr)
        labels = clusterer.labels_

        a = np.zeros(d_array.shape)
        for clust, tup in zip(labels, in_arr):
            if clust >= 0:
                a[tuple(tup)] = clust + 1
            else:
                a[tuple(tup)] = clust
        return a

