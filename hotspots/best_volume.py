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

The :mod:`hotspots.extraction` module contains classes to extract
valuable information from the calculated Fragment Hotspot Maps.

The main fields can be broken due using a fixed volume or fixed score mode.

This approach enables a ranking of ligand sized cavity descriptions to be
extracted. For example, these grids can be used to generate pharmacophores.

The main classes of the :mod:`hotspots.best_volume` module are:
    -HotspotFeatures
    -ExtractHotspot
    -HotspotCreator
"""
from __future__ import print_function
from os.path import join, dirname

import hotspot_calculation # import HotspotResults, _RunSuperstar
from grid_extension import Grid
from scipy import optimize
from skimage import feature
from hotspot_utilities import Helper


class HotspotResults(hotspot_calculation.HotspotResults):
    """
    A class to handle the construction of indivdual hotspots
    """
    class Optimiser(object):
        """
        A class to handle the optimisation operations
        """
        def __init__(self, mask, settings, peak=None):
            self.peak = peak
            self.mask = mask
            self.settings = settings

        def _count_island_points(self, threshold):
            """
            For a given island, the difference between the target number of grid points and the actual number of
             grid points is returned
            :param threshold:
            :return: int
            """
            island = self.mask.get_best_island(threshold, mode="score", peak=self.peak)
            if island is None:
                return 999999
            points = (island > threshold).count_grid()
            if self.peak:
                padding = self.settings.padding_factor * self.settings.num_gp
                return abs((self.settings.num_gp + padding) - points)
            else:
                return abs(self.settings.num_gp - points)

        def _count_grid_points(self, threshold):
            """
            For a given island, the difference between the target number of grid points and the actual number of
             grid points is returned
            :param threshold:
            :return: int
            """
            points = (self.top_island > threshold).count_grid()
            return abs(self.settings.num_gp - points)

        def _reselect_points(self, threshold):
            """
            looks within the top islands bounding box for other points above threshold.
            :param threshold: float, island threshold
            :return:
            """
            self.top_island = self.mask.get_best_island(threshold=threshold, mode="score", peak=self.peak)
            new_threshold = optimize.fminbound(self._count_grid_points, 0, 30, xtol=0.025)
            best_island = (self.top_island > new_threshold) * self.top_island

            return new_threshold, best_island

        def optimize_island_threshold(self):
            """
            Takes the input mask and finds the island threshold which returns the desired volume
            :return:
            """
            threshold = optimize.fminbound(self._count_island_points, 0, 30, xtol=0.025)
            if threshold > 29:
                threshold = 1
            best_island = self.mask.get_best_island(threshold=threshold, mode='score', peak=self.peak)
            best_island = (best_island > threshold) * best_island

            if best_island.count_grid() > self.settings.num_gp:
                threshold += 0.01
                best_island = (best_island > threshold) * best_island

            # new_threshold, best_island = self._reselect_points(threshold=threshold)
            print("target = {}, actual = {}".format(self.settings.num_gp, best_island.count_grid()))

            return threshold, best_island


    # def __init__(self, single_grid, mask_dic, settings, prot, large_cavities, superstar=None, seed=None):
    #     """
    #
    #     :param single_grid:
    #     :param mask_dic:
    #     :param settings:
    #     :param seed:
    #     """
    #     self._single_grid = single_grid
    #     self._mask_dic = mask_dic
    #     self.settings = settings
    #     self._protein = prot
    #
    #     blank = self.single_grid.copy_and_clear()
    #     if self.settings.mode == "seed":
    #         self._seed = seed
    #         self._seed_value = single_grid.value_at_point(self.seed)
    #         self.sphere = blank
    #         self.sphere.set_sphere(point=self.seed,
    #                                radius=settings.search_radius,
    #                                value=1,
    #                                scaling='None'
    #                                )
    #         self._mask = single_grid * self.sphere
    #
    #     else:
    #         self._mask = single_grid
    #
    #     self._rank = None
    #     self._threshold, self.best_island = self.optimize_island_threshold()
    #     self._location, self._features = self.get_interaction_type()
    #     self.rank_features()
    #     self._score = self.score_hotspot(mode="apolar_peaks")
    #     self.hotspot_result = self.get_hotspot_result()
    #     # if large_cavities:
    #     #     self.elaboration_potential = self.get_elaboration_potential(large_cavities=large_cavities)
    #     # else:
    #     #     self.elaboration_potential = None
    #
    #     if superstar:
    #         self.superstar_results = self.get_superstar_result(superstar)
    #         self.calc_feature_profile()

    # @staticmethod
    # def from_file(hr):
    #     location = hr.super_grids["apolar"].centroid()
    #     best_island = Grid.get_single_grid(hr.super_grids)
    #     # peak = best_island.centroid()
    #     # features
    #     #threshold
    #     #score
    #     #features =
    #
    #
    #     return ExtractedHotspot(hotspot_result=hr,
    #                             location=location,
    #                             best_island=best_island)

    @staticmethod
    def from_hotspot(single_grid, mask_dic, settings, prot, seed=None):
        """
        create a Extracted Hotspot object from HotspotResult object
        :param single_grid:
        :param mask_dic:
        :param settings:
        :param prot:
        :param seed:
        :return:
        """

        mask = single_grid
        optimiser = HotspotResults.Optimiser(mask=mask, settings=settings, peak=seed)
        threshold, best_island = optimiser.optimize_island_threshold()

        location, features = HotspotResults.get_interaction_type(mask_dic, best_island, threshold, settings)
        grd_dict = HotspotResults.get_grid_dict(location, features, settings)

        hr = HotspotResults(super_grids=grd_dict, protein=prot)
        hr.threshold = threshold
        hr.best_island = best_island
        hr.location = location
        hr.features = features
        hr.score = hr.score()
        hr.rank = hr._rank_features()
        return hr

    # def count_island_points(self, threshold):
    #     """
    #     For a given island, the difference between the target number of grid points and the actual number of
    #      grid points is returned
    #     :param threshold:
    #     :return: int
    #     """
    #     island = self.mask.get_best_island(threshold, mode="score")
    #     if island is None:
    #         return 999999
    #     points = (island > threshold).count_grid()
    #     # EXTRA PADDING FOR RUN
    #     padding = self.settings.padding_factor * self.settings.num_gp
    #
    #     return abs((self.settings.num_gp + padding) - points)
    #
    # def count_grid_points(self, threshold):
    #     """
    #     For a given island, the difference between the target number of grid points and the actual number of
    #      grid points is returned
    #     :param threshold:
    #     :return: int
    #     """
    #     points = (self.top_island > threshold).count_grid()
    #
    #     return abs(self.settings.num_gp - points)
    #
    # def reselect_points(self, threshold):
    #     """
    #     looks within the top islands bounding box for other points above threshold.
    #     :param threshold: float, island threshold
    #     :return:
    #     """
    #     self.top_island = self.mask.get_best_island(threshold=threshold, mode="score")
    #     # give some extra padding
    #     #step_out = self.top_island.step_out_mask(nsteps=2)
    #
    #     # self.top_island = Grid.super_grid(1, *[i for i in self.top_island.islands(threshold=threshold)])
    #
    #     new_threshold = optimize.fminbound(self.count_grid_points, 0, 30, xtol=0.025, disp=1)
    #
    #     best_island = (self.top_island > new_threshold)*self.top_island
    #
    #     return new_threshold, best_island
    #
    # def optimize_island_threshold(self):
    #     """
    #     Takes the input mask and finds the island threshold which returns the desired volume
    #     :return:
    #     """
    #     threshold = optimize.fminbound(self.count_island_points, 0, 30, xtol=0.025, disp=1)
    #     if threshold > 29:
    #         threshold = 1
    #
    #     if self.settings.mode == "seed":
    #         new_threshold, best_island = self.reselect_points(threshold=threshold)
    #         print "target = {}, actual = {}".format(self.settings.num_gp, best_island.count_grid())
    #         return new_threshold, best_island
    #
    #     else:
    #         new_threshold, best_island = self.reselect_points(threshold=threshold)
    #         print "target = {}, actual = {}".format(self.settings.num_gp, best_island.count_grid())
    #         return new_threshold, best_island

    @staticmethod
    def get_interaction_type(mask_dic, best_island, threshold, settings):
        """
        seperates single grid into grid by interaction type
        :return:
        """
        features = []
        location = None

        common_best_island = mask_dic["apolar"].common_boundaries(best_island)
        features_in_vol = {p: g * (g & common_best_island) for p, g in mask_dic.items()}
        location = features_in_vol["apolar"]
        features = HotspotResults._get_features(interaction_dict=features_in_vol,
                                                threshold=threshold,
                                                min_feature_gp=settings.min_feature_gp)

        return location, features

        #
        # for probe, g in mask_dic.items():
        #     if probe == "apolar":
        #         location = (g & common_best_island) * g
        #         if location.count_grid() == 0:
        #             raise RuntimeError("No apolar location found")
        #     else:
        #         features_in_vol = g * (g & common_best_island)
        #         if len(features_in_vol.islands(threshold=threshold)) > 0:
        #             features.extend(HotspotResults._get_features(probe, features_in_vol, threshold, settings))
        # return location, features

    # @staticmethod
    # def _get_features(probe, g, threshold, settings):
    #     """
    #     returns Hotspot Feature object with a score to enable ranking
    #     :param probe:
    #     :param g:
    #     :return:
    #     """
    #     return [HotspotResults.HotspotFeature(probe, island) for island in g.islands(threshold=threshold)
    #             if island.count_grid() > settings.min_feature_gp]

    # def rank_features(self):
    #     """
    #     rank features based upon feature score (TO DO: modify score if required)
    #     :return:
    #     """
    #     feature_by_score = {feat.score: feat for feat in self.features}
    #     score = sorted([f[0] for f in feature_by_score.items()], reverse=True)
    #     for i, key in enumerate(score):
    #         feature_by_score[key]._rank = int(i + 1)

    # def score_hotspot(self, mode="apolar_peaks"):
    #     """
    #     rank extracted hotspot. Modes:
    #         -apolar_peaks
    #         -all_peaks
    #     :return:
    #     """
    #     if self.settings.mode == "seed":
    #         apolar_score = self.seed_value
    #         self.alternative_score = np.mean([feat.score for feat in self.features] + [self.seed_value])
    #     else:
    #         islands = self.location.islands(threshold=self.threshold)
    #         i = [island.value_at_point(island.centroid()) for island in islands]
    #         apolar_score = sum(i)/len(i)
    #
    #     if mode == "apolar_peaks":
    #         return apolar_score
    #
    #     elif mode == "all_peaks":
    #             vals = [feat.score for feat in self.features] + [apolar_score]
    #             score = np.mean(vals)
    #     else:
    #         raise RuntimeError("Mode not support, request feature")
    #
    #     return score

    @staticmethod
    def get_grid_dict(location, features, settings):
        """
        Creates super grid dict from location and features
        :param location:
        :param features:
        :param settings:
        :return:
        """
        grid_dic = {"apolar": location.minimal()}
        interaction_types = set([feat.feature_type for feat in features])
        feature_by_score = {f.score: f for f in features}
        features = [feature_by_score[s]
                    for s in sorted([f[0] for f in feature_by_score.items()], reverse=True)][:settings.max_features - 1]
        for probe in interaction_types:
            if settings.mode == "seed":
                grids = [feat.grid for feat in features
                         if feat.feature_type == probe and
                         feat.score >= settings.cutoff]

            else:
                grids = [feat.grid for feat in features if feat.feature_type == probe]
            if len(grids) == 0:
                grids = [location.minimal().copy_and_clear()]

            grid_dic.update({probe: Grid.super_grid(1, *grids)})

        return grid_dic

    # def get_hotspot_result(self):
    #     """
    #     formats the location and set of interaction feature into super_grids, a dict of interaction type: grid
    #     :return:
    #     """
    #
    #     apolar_min = self.location.minimal()
    #     grid_dic = {"apolar": apolar_min}
    #     for probe in self.mask_dic.keys():
    #         if probe != "apolar":
    #             if self.settings.mode == "seed":
    #                 grids = [feat.grid for feat in self.features
    #                          if feat.feature_type == probe and
    #                          feat.rank <= self.settings.max_features and
    #                          feat.score >= self.settings.cutoff]
    #             else:
    #                 grids = [feat.grid for feat in self.features
    #                          if feat.feature_type == probe]
    #
    #             if len(grids) == 0:
    #                 grids = [apolar_min.copy_and_clear()]
    #
    #             grid_dic.update({probe: Grid.super_grid(1, *grids)})
    #         else:
    #             continue
    #
    #     return HotspotResults(grid_dict=grid_dic,
    #                           protein=self.prot,
    #                           sampled_probes=None,
    #                           buriedness=None)

        
    # def get_superstar_result(self, superstar_results):
    #     """
    #     finds the overlap between the extracted hotspot and the superstar results
    #     :param superstar_result:
    #     :return:
    #     """
    #     # TO DO: ALLOW SS RUN IN EXTRACTED_HOTSPOT CLASS
    #
    #     extracted_superstar = []
    #
    #     for result in superstar_results:
    #         common_best_island, common_result_grid = Grid.common_grid(self.best_island, result.grid)
    #         ss_boundary = (common_best_island & common_result_grid) * common_result_grid
    #         new = copy.copy(result)
    #         if len(ss_boundary.islands(threshold=2)) != 0:
    #             g = Grid.super_grid(2, *ss_boundary.islands(threshold=2))
    #             threshold = g.grid_score(threshold=1, percentile=50)
    #             print threshold
    #             new.grid = (g > threshold) * g
    #         else:
    #             new.grid = ss_boundary.copy_and_clear()
    #
    #         extracted_superstar.append(new)
    #
    #     return extracted_superstar
    #
    # def calc_feature_profile(self):
    #     """
    #     for each hotspot feature, the overlap between the feature sphere and superstar result is calculated
    #     this is stored as an HotspotFeature attribute (superstar profile)
    #     :return:
    #     """
    #     for feat in self.features:
    #         super_profile = []
    #
    #         for result in self.superstar_results:
    #             common_result_grid, common_sphere = Grid.common_grid(result.grid, feat.sphere)
    #             super_sphere = (common_sphere & common_result_grid) * common_result_grid
    #
    #             if len(super_sphere.islands(threshold=2)) != 0:
    #                 result.grid = Grid.super_grid(2, *super_sphere.islands(threshold=2))
    #
    #             else:
    #                 result.grid = feat.sphere.copy_and_clear()
    #
    #             super_profile.append(result)
    #
    #         feat.superstar_profile = super_profile


class Extractor(object):
    """
    A class to handle the extraction of discrete, fragment size hotspots from the original maps
    """

    def __init__(self, hr, settings=None, mode="seed", volume="125", pharmacophores=True):
        """

        :param hr:
        :param out_dir:
        :param settings:
        """
        print("Initialise Extractor")
        self.hotspot_result = hr
        if settings is None:
            self.settings = self.Settings()
        else:
            self.settings = settings

        self.settings.mode = mode
        self.settings.volume = volume
        self.settings.pharmacophore = pharmacophores
        self.out_dir = None

        # fragment hotspot post processing
        if self.settings.mode == "seed":
            hr.super_grids = self.grid_post_process(hr.super_grids)

        else:
            hr.super_grids.update({probe: g.max_value_of_neighbours() for probe, g in hr.super_grids.items()})

        try:
            hr.super_grids["negative"] = hr.super_grids["negative"].deduplicate(hr.super_grids["acceptor"],
                                                                                threshold=10,
                                                                                tolerance=2)

            hr.super_grids["positive"] = hr.super_grids["positive"].deduplicate(hr.super_grids["donor"],
                                                                                threshold=10,
                                                                                tolerance=2)
        except KeyError:
            pass

        # enable extraction to run in seeded or global modes
        if self.settings.mode == "seed":
            print("    Mode: 'seed'")
            self._peaks = self.get_peaks()

        elif self.settings.mode == "global":
            print("    Mode: 'global'")
            self._peaks = None

        else:
            raise IOError("Mode not currently supported")

        # runs and ranks extraction procedure
        print("Generate Single Grid")
        self._masked_dic, self._single_grid = Grid.get_single_grid(self.hotspot_result.super_grids)
        # self._large_cavities = self._get_large_cavities()
        print("Extracting Hotspots")
        self.extracted_hotspots = self._get_extracted_hotspots()

        # remove HS with > 80% overlap
        # for i, h in enumerate(self.extracted_hotspots):
        #     for i in range(i, len(self.extracted_hotspots) + 1):
        #         overlap = h.best_island.overlap(self.extracted_hotspots[i])
        #         if overlap > 0.8


        self._rank_extracted_hotspots()

        # generate pharmacophores
        if self.settings.pharmacophore:
            self.get_pharmacophores()


    class Settings(object):
        """
        Default settings for hotspot extraction
        """

        def __init__(self):
            self.volume = 150
            self.cutoff = 14
            self.search_radius = int(5)
            self.mode = "seed"
            self.padding_factor = 0

            self.spacing = 0.5

            self.min_feature_gp = 5
            self.max_features = 10
            self.min_distance = 6
            self.island_max_size = 100
            self.sigma = 0.3

            self.drug_volume = 300
            self.buriedness_value = 4
            self.fragments = None
            self.lead = None

            self.pharmacophore = True

        @property
        def num_gp(self):
            return int(float(self.volume) / self.spacing ** 3)

    @property
    def single_grid(self):
        return self._single_grid

    @property
    def masked_dic(self):
        return self._masked_dic

    @property
    def peaks(self):
        return self._peaks

    def grid_post_process(self, super_grids):
        """
        carry out post-processing of fragment hotspot maps

        Limit the size of polar islands. Keep top scores upto X grid points
        :return:
        """
        for probe, g in super_grids.items():
            if probe == "apolar":
                super_grids.update({probe: g.max_value_of_neighbours()})

            else:
                h = g.max_value_of_neighbours()
                h = h.limit_island_size(self.settings.island_max_size)
                if h.bounding_box != super_grids["apolar"].bounding_box:
                    h = super_grids["apolar"].common_boundaries(g)

                super_grids.update({probe: h})

        # try:
        #     super_grids["negative"] = super_grids["negative"].deduplicate(super_grids["acceptor"],
        #                                                                         threshold=10,
        #                                                                         tolerance=2)
        #
        #     super_grids["positive"] = super_grids["positive"].deduplicate(super_grids["donor"],
        #                                                                         threshold=10,
        #                                                                         tolerance=2)
        # except KeyError:
        #     pass

        return super_grids

    def get_peaks(self):
        """
        find peak coordinates in apolar maps, used as seeds to find top volumes
        :return:
        """
        apolar = self.hotspot_result.super_grids["apolar"]
        peaks = feature.peak_local_max(apolar.get_array(),
                                       min_distance=self.settings.min_distance,
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
                print(peak)
                e = HotspotResults.from_hotspot(self.single_grid,
                                                self.masked_dic,
                                                self.settings,
                                                self.hotspot_result.protein,
                                                seed=peak)

                print(e.threshold)
                if e.threshold > 12:
                    extracted_hotspots.append(e)

        else:
            e = HotspotResults.from_hotspot(self.single_grid,
                                            self.masked_dic,
                                            self.settings,
                                            self.hotspot_result.protein,
                                            seed=None)

            extracted_hotspots.append(e)

        return extracted_hotspots

    def _rank_extracted_hotspots(self):
        """
        assigns rank based upon extracted hotspot score
        :return:
        """
        hotspot_by_score = {hotspot.score: hotspot for hotspot in self.extracted_hotspots}
        score = sorted([f[0] for f in hotspot_by_score.items()], reverse=True)
        print(score)

        for i, key in enumerate(score):
            hotspot_by_score[key].rank = int(i + 1)

        extracted_hotspots_by_rank = {h.rank: h for h in self.extracted_hotspots}
        self.extracted_hotspots = [value for (key, value) in sorted(extracted_hotspots_by_rank.items())]

        for i, hs in enumerate(self.extracted_hotspots):
            hs.identifier = "rank_{}".format(hs.rank)
            print("rank", hs.rank, "score", hs.score)

    def _select_cavity_grids(self, cavs):
        """get empty cavity grids"""
        grds = [Grid(origin=cav.bounding_box[0],
                     far_corner=cav.bounding_box[1],
                     spacing=self.settings.spacing,
                     default=0)
                for cav in cavs]

        if self.settings.mode == "seed":
            filtered = set([g for seed in [p for p in self.peaks]
                            for g in grds
                            if g.contains_point(seed)])

        else:
            raise IOError("Currently only available in seed mode")

        return filtered

    def get_pharmacophores(self):
        """
        generates a pharmacophore model, stores as attribute of hotspot result
        :return:
        """
        for i, hotspot in enumerate(self.extracted_hotspots):
            hotspot.pharmacophore = hotspot.get_pharmacophore_model(identifier=hotspot.identifier)


    def _write(self, out_dir, mode="best_islands"):
        """
        write out information to aid debugging: valid modes:
            -peaks:
            -locations: spheres and islands at apolar peak locations
            -features: islands and probes at feature point locations
        """

        if mode == "peaks":
            out_dir = Helper.get_out_dir(join(out_dir))
            pymol_out = 'from pymol import cmd\nfrom pymol.cgo import *\n'
            for i, peak in enumerate(self.peaks):
                score = "{0:.2f}".format(self.hotspot_result.super_grids["apolar"].value_at_point(peak))
                sphere = 'score_{0} = [COLOR, 1.00, 1.000, 0.000] + ' \
                         '[ALPHA, 0.8] + ' \
                         '[SPHERE, float({1}), float({2}), float({3}), float(0.5)]\n' \
                    .format(i, peak[0], peak[1], peak[2])
                pymol_out += sphere
                pymol_out += '\ncmd.load_cgo(score_{1}, "score_{0}", 1)'.format(score, i)
                pymol_out += '\ncmd.group("Peaks", members="score_{0}")\n'.format(score)
            with open(join(out_dir, "peaks.py"), "w") as pymol_file:
                pymol_file.write(pymol_out)

        elif mode == "best_islands":
            out_dir = Helper.get_out_dir(join(out_dir, "best_islands"))
            pymol_out = 'from pymol import cmd\nfrom pymol.cgo import *\n'
            thresholds = []
            for i, extracted in enumerate(self.extracted_hotspots):
                extracted.best_island.write(join(out_dir, "island_{}.grd".format(i)))
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

            with open(join(dirname(out_dir), "best_islands.py"), "w") as pymol_file:
                pymol_file.write(pymol_out)

        else:
            raise IOError("mode not supported")


