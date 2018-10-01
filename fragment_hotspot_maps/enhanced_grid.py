'''
The main class of the :mod:`ccdc_development.Grid`.

This is an internal extension of :class:`ccdc.grid.Grid` that adds in potential new features
for review and development. For example, gaussian smoothing function:

.. code-block:: python

>>> from ccdc_development import Grid
>>> from scipy import ndimage
>>> import numpy as np
>>>
>>> grd = Grid.from_file(<path/to/file>)
>>> new_grd = grd.gaussian(sigma=0.3)
>>>

I.e. this returns a new grid object in which gaussian smoothing has been applied
'''
#############################################################################
from ccdc.utilities import Grid
import numpy as np
import numpy.ma as ma
from scipy import ndimage


def grid_score(self, threshold=0, percentile=75):
    """
    take a group and return average score of points above threshold
    :param g:
    :param threshold
    :param percentile

    :return:
    """
    array = self.get_array()
    masked_array = ma.masked_less_equal(array, threshold)
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
    grid = self.copy()
    grid *= 0
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
    if self.value_at_point(point) > threshold:
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
    grid = self.copy()
    grid *= 0
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

    return Grid.super_grid(2, *grid.islands(threshold=1))


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
    all_islands = set([jsland for jsland in self.islands(threshold=threshold)])
    bin_islands = set([jsland for jsland in all_islands
                       for island in major.islands(threshold=threshold)
                       if jsland.contains_point(island.centroid(), tolerance=tolerance)])

    retained_jslands = list(all_islands - bin_islands)
    print retained_jslands

    if len(retained_jslands) == 0:
        blank = self.copy() * 0
        return blank
    else:
        return Grid.super_grid(1, *retained_jslands)


def copy_and_clear(self):
    """
    make a new empty grid
    :param self:
    :return:
    """

    g = self.copy()
    g *= 0
    return g

# def _common_grid(self, g1, g2, padding=1):
#     """
#     finds a common frame of reference for two grids
#
#     :param g1: a :class: `ccdc.utilities.Grid` instance
#     :param g2: a :class: `ccdc.utilities.Grid` instance
#     :param padding: int, additional grid point in x, y, z directions
#     :return:
#     """
#
#     sg = Grid.super_grid(padding, g1, g2)
#     out_g = self._copy_and_clear(sg)
#     out1 = Grid.super_grid(padding, g1, out_g)
#     out2 = Grid.super_grid(padding, g2, out_g)
#     return out1, out2

Grid.grid_score = grid_score
Grid.indices_to_point = indices_to_point
Grid.point_to_indices = point_to_indices
Grid.gaussian = gaussian
Grid.contains_point = contains_point
Grid.get_array = get_array
Grid.restricted_volume = restricted_volume
Grid.centroid = centroid
Grid.deduplicate = deduplicate
Grid.copy_and_clear = copy_and_clear
