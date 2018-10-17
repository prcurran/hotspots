#
# This code is Copyright (C) 2015 The Cambridge Crystallographic Data Centre
# (CCDC) of 12 Union Road, Cambridge CB2 1EZ, UK and a proprietary work of CCDC.
# This code may not be used, reproduced, translated, modified, disassembled or
# copied, except in accordance with a valid licence agreement with CCDC and may
# not be disclosed or redistributed in any form, either in whole or in part, to
# any third party. All copies of this code made in accordance with a valid
# licence agreement as referred to above must contain this copyright notice.
#
# No representations, warranties, or liabilities are expressed or implied in the
# supply of this code by CCDC, its servants or agents, except where such
# exclusion or limitation is prohibited, void or unenforceable under governing
# law.
#
"""

The :mod:`fragment_hotspot_maps.extraction` module contains classes to extract
valuable information from the calculated Fragment Hotspot Maps.

The main fields can be broken due using a fixed volume or fixed score mode.

This approach enables a ranking of ligand sized cavity descriptions to be
extracted. For example, these grids can be used to generate pharmacophores.

The main classes of the :mod:`fragment_hotspot_maps.best_volume` module are:
    -HotspotFeatures
    -ExtractHotspot
    -HotspotCreator
"""
from enhanced_grid import Grid
from os.path import join, exists
import operator
from skimage import feature
from scipy import optimize
import numpy as np
from os import mkdir

from pprint import pprint


class HotspotFeatures(object):
    """
    class to hold polar islands above threshold "features"
    purpose: enables feature ranking
    """
    def __init__(self, feature_type, grid):
        self._feature_type = feature_type
        self._grid = grid

        self._coordinates = grid.centroid()
        self._count = grid.count_grid()

        self._score = self.score_feature()
        self._rank = None

    @property
    def feature_type(self):
        return self._feature_type

    @property
    def grid(self):
        return self._grid

    @property
    def coordinate(self):
        return self._coordinates

    @property
    def count(self):
        return self._count

    @property
    def score(self):
        return self._score

    @property
    def rank(self):
        return self._rank

    def score_feature(self):
        nx, ny, nz = self.grid.nsteps

        return sum([self.grid.value(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)])/self._count


class ExtractedHotspot(object):
    """
    A class to handle the construction of indivdual hotspots
    """

    def __init__(self, single_grid, mask_dic, settings, seed=None):
        """

        :param single_grid:
        :param mask_dic:
        :param settings:
        :param seed:
        """

        self._single_grid = single_grid
        self._mask_dic = mask_dic
        self.settings = settings
        blank = self.single_grid.copy_and_clear()
        print self.settings.mode
        if self.settings.mode == "seed":
            self._seed = seed
            self._seed_value = single_grid.value_at_point(self.seed)
            self.sphere = blank
            self.sphere.set_sphere(point=self.seed,
                                   radius=settings.search_radius,
                                   value=1,
                                   scaling='None'
                                   )
            self._mask = single_grid * self.sphere

        else:
            self._mask = single_grid

        self.super_grids = None
        self._rank = None

        self._threshold, self.best_island = self.optimize_island_threshold()
        self._location, self._features = self.get_interaction_type()
        self.rank_features()
        self._score = self.score_hotspot(mode="apolar_peaks")

    @property
    def seed(self):
        return self._seed

    @property
    def seed_value(self):
        return self._seed_value

    @property
    def mask(self):
        return self._mask

    @property
    def single_grid(self):
        return self._single_grid

    @property
    def mask_dic(self):
        return self._mask_dic

    @property
    def threshold(self):
        return self._threshold

    @property
    def location(self):
        return self._location

    @property
    def features(self):
        return self._features

    @property
    def rank(self):
        return self._rank

    @property
    def score(self):
        return self._score

    def count_island_points(self, threshold):
        """
        For a given island, the difference between the target number of grid points and the actual number of
         grid points is returned
        :param threshold:
        :return: int
        """
        island = self.mask.get_best_island(threshold, mode="score")
        if island is None:
            return 999999
        points = (island > threshold).count_grid()

        return abs(self.settings.num_gp - points)

    def count_grid_points(self, threshold):
        """
        For a given island, the difference between the target number of grid points and the actual number of
         grid points is returned
        :param threshold:
        :return: int
        """
        points = (self.top_island > threshold).count_grid()

        return abs(self.settings.num_gp - points)

    def reselect_points(self, threshold):
        """
        looks within the top islands bounding box for other points above threshold.
        :param threshold: float, island threshold
        :return:
        """
        self.top_island = self.mask.get_best_island(threshold=threshold, mode="score")
        print self.top_island.count_grid()
        new_threshold = optimize.fminbound(self.count_grid_points, 0, 30, xtol=0.025, disp=3)
        self.top_island.write("C:/Users/pcurran/Desktop/test_dev_inputs/hotspot_volume/top_island.grd")
        best_island = self.top_island > new_threshold
        print best_island.count_grid()
        return new_threshold, best_island,

    def optimize_island_threshold(self):
        """
        Takes the input mask and finds the island threshold which returns the desired volume
        :return:
        """

        threshold = optimize.fminbound(self.count_island_points, 0, 30, xtol=0.025, disp=3)
        print threshold
        if threshold > 29:
            threshold = 1

        if self.settings.mode == "seed":
            best_island = self.mask.get_best_island(threshold=threshold, mode="score")
            return threshold, best_island
        else:
            new_threshold, best_island = self.reselect_points(threshold=threshold)
            print new_threshold
            g = best_island > new_threshold
            print "target = {}, actual = {}".format(self.settings.num_gp, g.count_grid())
            return new_threshold, best_island

    def get_interaction_type(self):
        """
        seperates single grid into grid by interaction type
        :return:
        """
        features = []
        location = None

        common_best_island = self.mask_dic["apolar"].common_boundaries(self.best_island)

        for probe, g in self.mask_dic.items():
            if probe == "apolar":
                location = (g & common_best_island) * g

            else:
                features_in_vol = g * (g & common_best_island)
                if len(features_in_vol.islands(threshold=self.threshold)) > 0:
                    features.extend(self._get_features(probe, features_in_vol))
        return location, features

    def _get_features(self, probe, g):
        """
        returns Hotspot Feature object with a score to enable ranking
        :param probe:
        :param g:
        :return:
        """
        return [HotspotFeatures(probe, island) for island in g.islands(threshold=self.threshold)
                if island.count_grid() > self.settings.min_feature_gp]

    def rank_features(self):
        """
        rank features based upon feature score (TO DO: modify score if required)
        :return:
        """
        feature_by_score = {feat.score: feat for feat in self.features}
        score = sorted([f[0] for f in feature_by_score.items()], reverse=True)
        for i, key in enumerate(score):
            feature_by_score[key]._rank = int(i + 1)

    def score_hotspot(self, mode="apolar_peaks"):
        """
        rank extracted hotspot. Modes:
            -apolar_peaks
            -all_peaks
        :return:
        """
        if self.settings.mode == "seed":
            apolar_score = self.seed_value
            self.alternative_score = np.mean([feat.score for feat in self.features] + [self.seed_value])
        else:
            islands = self.location.islands(threshold=self.threshold)
            i = [island.value_at_point(island.centroid()) for island in islands]
            apolar_score = sum(i)/len(i)

        if mode == "apolar_peaks":
            return apolar_score

        elif mode == "all_peaks":
                vals = [feat.score for feat in self.features] + [apolar_score]
                score = np.mean(vals)
        else:
            raise RuntimeError("Mode not support, request feature")

        return score


class HotspotCreator(object):
    """
    A class to handle the extraction of discrete, fragment size hotspots from the original maps
    """

    def __init__(self, super_grids, out_dir, settings):
        import time
        self.start = time.time()
        self.settings = settings
        self.out_dir = out_dir
        self.super_grids = {}

        for probe, g in super_grids.items():

            if g.bounding_box != super_grids["apolar"].bounding_box:
                g = super_grids["apolar"].common_boundaries(g)

            self.super_grids.update({probe: g.gaussian(0.25).max_value_of_neighbours()})

        try:
            self.super_grids["negative"] = self.super_grids["negative"].deduplicate(self.super_grids["acceptor"],
                                                                                    threshold=14,
                                                                                    tolerance=3)

            self.super_grids["positive"] = self.super_grids["positive"].deduplicate(self.super_grids["donor"],
                                                                                    threshold=14,
                                                                                    tolerance=3)
        except KeyError:
            pass

        if self.settings.mode == "seed":
            self._peaks = self.get_peaks()
            self.write(mode="peaks")

        elif self.settings.mode == "global":
            self._peaks = None

        else:
            raise IOError("Mode not currently supported")

        self._masked_dic, self._single_grid = self.get_single_grid()
        self.extracted_hotspots = self._get_extracted_hotspots()
        self._rank_extracted_hotspots()
        self._format_data()

        self.write(mode="best_islands")

    class Settings(object):
        """
        Default settings for hotspot extraction
        """

        def __init__(self):
            self.volume = 150
            self.spacing = 0.5
            self.num_gp = int(float(self.volume) / self.spacing**3)
            self.cutoff = 14
            self.search_radius = int(5)
            self.mode = "seed"
            self.min_feature_gp = 5
            self.max_features = 10

    @property
    def single_grid(self):
        return self._single_grid

    @property
    def masked_dic(self):
        return self._masked_dic

    @property
    def peaks(self):
        return self._peaks

    def get_peaks(self):
        """
        find peak coordinates in apolar maps, used as seeds to find top volumes
        :return:
        """
        apolar = self.super_grids["apolar"]
        peaks = feature.peak_local_max(apolar.get_array(),
                                       min_distance=6,
                                       num_peaks_per_label=1,
                                       threshold_abs=self.settings.cutoff)
        peak_by_value = {}
        for peak in peaks:
            val = apolar.value(int(peak[0]), int(peak[1]), int(peak[2]))
            if val > self.settings.cutoff:
                if val in peak_by_value:
                    peak_by_value[val].append((peak[0], peak[1], peak[2]))
                else:
                    peak_by_value.update({val:[(peak[0], peak[1], peak[2])]})

        average_peaks = []
        for key in peak_by_value.keys():
            x = [point[0] for point in peak_by_value[key]]
            y = [point[1] for point in peak_by_value[key]]
            z = [point[2] for point in peak_by_value[key]]
            average_peaks.append(apolar.indices_to_point(int(sum(x)/len(x)),
                                                         int(sum(y)/len(y)),
                                                         int(sum(z)/len(z))

                                                         )
                                 )
        return average_peaks

    def _get_extracted_hotspots(self):
        """
        locate peaks in apolar maps and define fragment size volume
        :return: list of peak coordinates
        """
        extracted_hotspots = []
        if self.settings.mode == "seed":
            for peak in self.peaks:
                e = ExtractedHotspot(self.single_grid, self.masked_dic, self.settings, seed=peak)
                if e.threshold > 12:
                    extracted_hotspots.append(e)

        else:
            e = ExtractedHotspot(self.single_grid, self.masked_dic, self.settings, seed=None)
            extracted_hotspots.append(e)

        return extracted_hotspots

    def get_single_grid(self):
        """
        from a collection of grids, create one grid with the maximum value at each grid point
        :return: dictionary of mask by interaction type, single maximal grid
        """
        mask_dic = {}
        for probe, grid in self.super_grids.items():
            other_grids = [self.super_grids[p] for p in self.super_grids.keys() if p != probe]
            mask_dic.update({probe: grid * grid.multi_max_mask(other_grids)})
        blank = self.super_grids["apolar"].copy_and_clear()

        return mask_dic, reduce(operator.add, mask_dic.values(), blank)

    def _format_data(self):
        """
        create grid_dic and probe_dic to create hotspot results object
        """
        for extracted in self.extracted_hotspots:
            apolar_min = extracted.location.minimal()
            grid_dic = {"apolar": apolar_min}
            for probe in self.super_grids.keys():
                if probe != "apolar":
                    grids = [feat.grid for feat in extracted.features
                             if feat.feature_type == probe and
                             feat.rank <= self.settings.max_features and
                             feat.score >= self.settings.cutoff]
                    if len(grids) == 0:
                        grids = [apolar_min.copy_and_clear()]

                    grid_dic.update({probe: Grid.super_grid(1, *grids)})
                else:
                    continue
            extracted.super_grids = grid_dic

    def _rank_extracted_hotspots(self):
        """
        assigns rank based upon extracted hotspot score
        :return:
        """
        hotspot_by_score = {hotspot.score: hotspot for hotspot in self.extracted_hotspots}
        score = sorted([f[0] for f in hotspot_by_score.items()], reverse=True)
        print score
        for i, key in enumerate(score):
            hotspot_by_score[key]._rank = int(i + 1)

        extracted_hotspots_by_rank = {h.rank: h for h in self.extracted_hotspots}
        self.extracted_hotspots = [value for (key, value) in sorted(extracted_hotspots_by_rank.items())]
        for hs in self.extracted_hotspots:
            print "rank", hs.rank, "score", hs.score

    def write(self, mode='peaks'):
        """
        write out information to aid debugging: valid modes:
            -peaks:
            -locations: spheres and islands at apolar peak locations
            -features: islands and probes at feature point locations
        """
        pymol_out = 'from pymol import cmd\nfrom pymol.cgo import *\n'

        if mode == "peaks":
            for i, peak in enumerate(self.peaks):
                score = "{0:.2f}".format(self.super_grids["apolar"].value_at_point(peak))
                sphere = 'score_{0} = [COLOR, 1.00, 1.000, 0.000] + ' \
                         '[ALPHA, 0.8] + ' \
                         '[SPHERE, float({1}), float({2}), float({3}), float(0.5)]\n' \
                    .format(i, peak[0], peak[1], peak[2])
                pymol_out += sphere
                pymol_out += '\ncmd.load_cgo(score_{1}, "score_{0}", 1)'.format(score, i)
                pymol_out += '\ncmd.group("Peaks", members="score_{0}")\n'.format(score)

        elif mode == "best_islands":
            od = join(self.out_dir, "hotspot_volume", "best_islands")
            if not exists(od):
                mkdir(od)
            thresholds = []
            for i, extracted in enumerate(self.extracted_hotspots):
                extracted.best_island.write(join(od, "island_{}.grd".format(i)))
                thresholds.append(extracted.threshold)
            pymol_out += """
nh = {0}
thresholds = {1}
for n in range(nh):
    cmd.load(r'best_islands/island_%s.grd' % (n), 'apolar_%s' % (n))
    cmd.isosurface('surface_apolar_%s' % (n), 'apolar_%s' % (n), thresholds[n])
    cmd.set('transparency', 0.7, 'surface_apolar_%s' % (n))
    cmd.color('yellow', 'surface_apolar_%s' % (n))
for n in range(nh):
    cmd.group('hotspot_%s'%(n), members= 'surface_apolar_%s'%(n))
    cmd.group('hotspot_%s'%(n), members= 'apolar_%s'%(n))""" \
            .format(len(self.extracted_hotspots), thresholds)

        else:
            raise IOError("write mode not supported")

        with open(join(self.out_dir, "hotspot_volume", "{}.py".format(mode)), "w") as pymol_file:
            pymol_file.write(pymol_out)
