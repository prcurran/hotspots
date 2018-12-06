from __future__ import print_function, division

from grid_extension import Grid
import os
from os.path import join, dirname, basename
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
from scipy import ndimage
from fragment_hotspot_maps.fragment_hotspot_maps import Hotspots
import pickle


class GridEnsemble(object):
    """Class that handles a grid of tuples for Hotspots across multiple structures"""

    def __init__(self):
        # self.stem = None
        self.prot_name = None
        self.probe = None
        self.path_list = None
        self.tup_max_length = None
        self.results_array = None
        self.common_grid_origin = None
        self.common_grid_far_corner = None
        self.common_grid_nsteps = None
        self.out_dir = None
        self.spacing = None

    def common_grids_from_paths(self):
        """
        Gets the coordinates of a grid that can fit all other grids (to use as basis for self.results_array)
        :return: 
        """
        print('Making array grid')
        grid_list = [Grid.from_file(f) for f in self.path_list if self.probe in basename(f)]
        common_grids = Grid.common_grid(grid_list)
        assert (len(set([c.bounding_box for c in common_grids])) == 1), "Common grids don't have same frames"
        self.spacing = common_grids[0].spacing
        assert (common_grids[0].spacing == 0.5), "Grid spacing not 0.5"
        self.tup_max_length = len(grid_list)
        self.common_grid_origin = common_grids[0].bounding_box[0]
        self.common_grid_far_corner = common_grids[0].bounding_box[1]
        self.common_grid_nsteps = common_grids[0].nsteps

        return common_grids

    def common_grids_from_grid_list(self, grid_list):
        common_grids = Grid.common_grid(grid_list)
        assert (len(set([c.bounding_box for c in common_grids])) == 1), "Common grids don't have same frames"
        self.spacing = common_grids[0].spacing
        assert (common_grids[0].spacing == 0.5), "Grid spacing not 0.5"
        self.tup_max_length = len(grid_list)
        self.common_grid_origin = common_grids[0].bounding_box[0]
        self.common_grid_far_corner = common_grids[0].bounding_box[1]
        self.common_grid_nsteps = common_grids[0].nsteps

    def get_results_array(self, common_grids):
        """
        constructs the numpy array containing the list of tuples.
        Adjusts the coordinates based on the origins of the grids
        self.results_array = numpy array
        :return:
        """
        print('Making results array')

        results_array = np.zeros(self.common_grid_nsteps, dtype=tuple)
        rec_spacing = 1 / self.spacing

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

    def pickle_GridEnsemble(self):
        """
        Saves MegaGrids as pickles.
        :return: 
        """
        pickle.dump(self, open(join(self.out_dir, '{}_{}_GridEnsemble.p'.format(self.prot_name, self.probe)), 'wb'))

    def make_grid(self, values):
        """
        Makes a grid to store output of ranges, max, means, etc 
        :param values: 
        :return: 
        """
        grid = Grid(origin=self.common_grid_origin,
                    far_corner=self.common_grid_far_corner,
                    spacing=0.5,
                    default=0.0,
                    _grid=None)

        nonz = self.results_array.nonzero()
        as_triads = zip(*nonz)
        for (i, j, k), v in zip(as_triads, values):
            grid._grid.set_value(int(i), int(j), int(k), v)

        return grid

    def get_gridpoint_histograms(self):
        """
        Makes and saves histograms for each point in the results array
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

    def get_gridpoint_means(self):
        """
        For each point in the MegaGrid, calculates the difference in score between each point in the tuple and the mean of the tuple. 
        :return: Python list
        """
        ind_array = np.indices(self.results_array.shape)
        means = []

        def get_means(x, y, z):
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Number of zeros', num_zeros, 'out of ', self.tup_max_length)
                hist_arr = np.array(self.results_array[x][y][z])
                means.extend(list(hist_arr - np.mean(hist_arr)))

        vget_means = np.vectorize(get_means, otypes=[list])
        vget_means(ind_array[0], ind_array[1], ind_array[2])
        return means

    def get_gridpoint_max(self):
        """
        For each point in the MegaGrid, calculates the max of the tuple
        :return: 
        """
        ind_array = np.indices(self.results_array.shape)
        maxes = []

        def get_max(x, y, z):
            """
            Would be funnier if I knew a Max.
            """
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Number of zeros: ', num_zeros)
                hist_arr = np.array(self.results_array[x][y][z])
                maxes.append(max(hist_arr))

        vget_max = np.vectorize(get_max, otypes=[list])
        vget_max(ind_array[0], ind_array[1], ind_array[2])
        return maxes

    def plot_gridpoint_spread(self, means):
        '''
        For each point in the 3D grid, plots the difference in score between each point in the tuple and the mean of the tuple. 
        
        '''
        mean_arr = np.array(means)
        (mu, sigma) = norm.fit(mean_arr)
        n, bins, patches = plt.hist(means, bins=40, normed=1)
        # print(bins)
        # y_fit = np.random.normal(mu, sigma, np.shape(mean_arr))
        y = mlab.normpdf(bins, mu, sigma)
        a = plt.plot(bins, y, 'r--', linewidth=2)
        bins = 0.5 * (bins[1:] + bins[:-1])
        y_fit = mlab.normpdf(bins, mu, sigma)
        ss_res = np.sum((n - y_fit) ** 2)
        ss_tot = np.sum((n - np.mean(n)) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        plt.title('Mu: {}, Sigma: {}, R^2 : {} '.format(round(mu, 2), round(sigma, 2), round(r2, 4)))
        hist_name = self.prot_name + '_{}_gridpoint_spread_fit'.format(self.probe)
        plt.savefig(join(self.out_dir, hist_name))
        # plt.show()
        plt.close()

    def get_gridpoint_ranges(self):
        ind_array = np.indices(self.results_array.shape)
        ranges = []

        def get_ranges(x, y, z):
            if isinstance(self.results_array[x][y][z], tuple):
                num_zeros = self.tup_max_length - len(self.results_array[x][y][z])
                if num_zeros != 0:
                    print('Number of zeros: ', num_zeros)
                hist_arr = np.array(self.results_array[x][y][z])
                ranges.append(max(hist_arr) - min(hist_arr))

        vget_ranges = np.vectorize(get_ranges, otypes=[list])
        vget_ranges(ind_array[0], ind_array[1], ind_array[2])
        return ranges

    def plot_gridpoint_ranges(self, ranges):
        '''
        Plots the range of the tuple values for each point in the 3D grid
        :return: ranges = list of the tuple ranges at each point
        '''

        plt.hist(ranges, bins=40, normed=0)
        plt.title('Score ranges for {} {}'.format(self.prot_name, self.probe))
        hist_name = self.prot_name + '_{}_score_ranges'.format(self.probe)
        plt.savefig(join(self.out_dir, hist_name))
        # plt.show()
        plt.close()

    def make_frequency_grid(self):
        '''
        Makes a grid that stores how many times each point has been sampled
        :return: ccdc.utilities Grid object
        '''
        r_a_grid = Grid(origin=self.common_grid_origin, far_corner=self.common_grid_far_corner, spacing=self.spacing)

        for x in range(self.results_array.shape[0]):
            for y in range(self.results_array.shape[1]):
                for z in range(self.results_array.shape[2]):
                    if isinstance(self.results_array[x][y][z], tuple):
                        r_a_grid.set_value(x, y, z, len(self.results_array[x][y][z]))

        r_a_grid.write(join(self.out_dir, '{}_Ghecom_frequency_{}.ccp4'.format(self.prot_name, self.probe)))
        return r_a_grid

    def make_ranges_grid(self):
        '''
        Makes a grid that stores the ranges of values of each point
        :return: ccdc.utilities Grid object
        '''
        r_a_grid = Grid(origin=self.common_grid_origin, far_corner=self.common_grid_far_corner, spacing=self.spacing)

        for x in range(self.results_array.shape[0]):
            for y in range(self.results_array.shape[1]):
                for z in range(self.results_array.shape[2]):
                    if isinstance(self.results_array[x][y][z], tuple):
                        hist_arr = np.array(self.results_array[x][y][z])
                        r_a_grid.set_value(x, y, z, (max(hist_arr) - min(hist_arr)))

        r_a_grid.write(join(self.out_dir, '{}_Ghecom_ranges_{}.ccp4'.format(self.prot_name, self.probe)))
        return r_a_grid

    def make_max_grid(self):
        '''
        Makes a grid that stores the maximum of values of each point
        :return: ccdc.utilities Grid object
        '''
        r_a_grid = Grid(origin=self.common_grid_origin, far_corner=self.common_grid_far_corner, spacing=self.spacing)

        for x in range(self.results_array.shape[0]):
            for y in range(self.results_array.shape[1]):
                for z in range(self.results_array.shape[2]):
                    if isinstance(self.results_array[x][y][z], tuple):
                        hist_arr = np.array(self.results_array[x][y][z])
                        r_a_grid.set_value(x, y, z, (max(hist_arr)))

        r_a_grid.write(join(self.out_dir, '{}_max_{}.ccp4'.format(self.prot_name, self.probe)))
        return r_a_grid

    def make_mean_grid(self):
        '''
        Makes a grid that stores the mean of the sampled values at each point
        :return: ccdc.utilities Grid object
        '''
        r_a_grid = Grid(origin=self.common_grid_origin, far_corner=self.common_grid_far_corner, spacing=self.spacing)

        for x in range(self.results_array.shape[0]):
            for y in range(self.results_array.shape[1]):
                for z in range(self.results_array.shape[2]):
                    if isinstance(self.results_array[x][y][z], tuple):
                        t_a = np.array(self.results_array[x][y][z])
                        r_a_grid.set_value(x, y, z, np.mean(t_a))

        r_a_grid.write(join(self.out_dir, '{}_Ghecom_mean_{}.ccp4'.format(self.prot_name, self.probe)))
        return r_a_grid

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

    def from_hotspot_maps(self, path_list, out_dir, prot_name, probe_name):
        """
        Creates a MegaGrid from a number of Hotspot maps for a certain probe
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

        common_grids = self.common_grids_from_paths()
        self.get_results_array(common_grids)

        return self.make_max_grid()

    def from_grid_list(self, grid_list, out_dir, prot_name, probe_name):
        """
        
        :param grid_list: 
        :param out_dir: 
        :param prot_name: 
        :param probe_name: 
        :return: 
        """
        self.out_dir = out_dir
        self.prot_name = prot_name
        self.probe = probe_name
        common_grids = self.common_grids_from_grid_list(grid_list)
        self.get_results_array(common_grids)

        return self.make_max_grid()

    ####### No longer in use ##########
    class EnsembleToHotspot(object):

        """Class that integrates ensembles with the rest of the hotspots"""

        def __init__(self):
            self.prot_name = None
            # self.protein = None
            self.path_list = None
            self.out_dir = None
            self.probe_list = None

        def result_from_ensembles(self, path_list, prot_name, out_dir):
            """
            Makes a hotspot results object from 
            :param path_list: list, paths to aligned hotspot maps of the same protein
            :param prot_name: str, protein name 
            :param out_dir:  str, path to where output is stored
            :param charged: 
            :return: HotspotResults
            """
            probe_list = ["acceptor", "apolar", "donor", "positive", "negative"]

            grid_dic = {}

            for p in probe_list:
                paths = [path for path in path_list if p in path]
                print(len(paths))
                if len(paths) == 0:
                    print("No maps for {} probe".format(p))
                    continue
                ens = GridEnsemble()
                grid_dic[p] = ens.from_hotspot_maps(paths, out_dir, prot_name, p)

            return Hotspots.HotspotResults(grid_dic, protein=None, fname=None, buriedness=None, sampled_probes=None)

        def ensemble_from_results(self, res_list, prot_name, out_dir, charged=False):
            """

            :param res_list: list of Hotspot.Results
            :param prot_name: str
            :param out_dir: str
            :return: HotspotResults
            """
            if charged:
                probe_list = ["acceptor", "apolar", "donor", "positive", "negative"]
            else:
                probe_list = ["acceptor", "apolar", "donor", ]

            grid_dic = {}

            for p in probe_list:
                grid_list_p = [r.super_grids[p] for r in res_list]
                ens = GridEnsemble()
                grid_dic[p] = ens.from_grid_list(grid_list_p, out_dir, prot_name, p)


if __name__ == "__main__":
    """
    from glob import glob

    paths = glob(join("/home/jin76872/Desktop/Mih/Data/BAZ2B_BRD1_Oct/BRD1_max_maps/*.ccp4"))
    out = "/home/jin76872/Desktop/Mih/Data/18_11_19_voronoi"
    pname = "BRD1"
    print("doinf stuff")
    tohot = EnsembleToHotspot()
    hr = tohot.result_from_ensembles(paths, pname, out)
    """
