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

The :mod:`hotspots.best_volume` module contains classes to extract
valuable information from the calculated Fragment Hotspot Maps.

The main fields can be broken due using a fixed volume or fixed score mode.

This approach enables a ranking of ligand sized cavity descriptions to be
extracted. For example, these grids can be used to generate pharmacophores.

The main classes of the :mod:`hotspots.best_volume` module are:
    -Extractor
    """
from __future__ import print_function, division
from os.path import join, dirname

import calculation  # import HotspotResults, _RunSuperstar
from grid_extension import Grid
from scipy import optimize
from skimage import feature
from hs_utilities import Helper


class _Results(calculation.Results):
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
            return abs(self.settings._num_gp - points)
            # if self.peak:
            #     padding = self.settings.padding_factor * self.settings._num_gp
            #     return abs((self.settings._num_gp + padding) - points)
            # else:
            #     return abs(self.settings._num_gp - points)

        def _count_grid_points(self, threshold):
            """
            For a given island, the difference between the target number of grid points and the actual number of
             grid points is returned
            :param threshold:
            :return: int
            """
            points = (self.top_island > threshold).count_grid()
            return abs(self.settings._num_gp - points)

        def _reselect_points(self, threshold):
            """
            looks within the top islands bounding box for other points above threshold.
            :param threshold: float, island threshold
            :return:
            """
            self.top_island = self.mask.get_best_island(threshold=threshold, mode="score", peak=self.peak)
            new_threshold = optimize.fminbound(self._count_grid_points, 0, 30, xtol=0.01)
            best_island = (self.top_island > new_threshold) * self.top_island

            return new_threshold, best_island

        def optimize_island_threshold(self):
            """
            Takes the input mask and finds the island threshold which returns the desired volume
            :return:
            """
            threshold = optimize.fminbound(self._count_island_points, 0, 30, xtol=0.025)
            if threshold >28:
                threshold = 1
            best_island = self.mask.get_best_island(threshold=threshold, mode='score', peak=self.peak)

            # If threshold is close to zero, keep all grid points
            try:
                best_island = (best_island > threshold) * best_island
            except TypeError:
                best_island = self.mask

            if best_island.count_grid() > self.settings._num_gp:
                threshold += 0.01
                best_island = (best_island > threshold) * best_island

            # new_threshold, best_island = self._reselect_points(threshold=threshold)
            # if (self.settings._num_gp * 0.85) < best_island.count_grid() < (self.settings._num_gp * 1.15):
            print("target = {}, actual = {}".format(self.settings._num_gp, best_island.count_grid()))
            return threshold, best_island
            # else:
            #     return None, None

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
        if seed:
            sphere = single_grid.copy_and_clear()
            sphere.set_sphere(point=seed, radius=settings._search_radius, value=1, scaling='None')
            mask = (sphere & single_grid) * single_grid
        else:
            mask = single_grid

        optimiser = _Results.Optimiser(mask=mask, settings=settings, peak=seed)
        threshold, best_island = optimiser.optimize_island_threshold()

        print("oooh",threshold)
        if best_island is not None:
            location, features = _Results.get_interaction_type(mask_dic, best_island, threshold, settings)
            grd_dict = _Results.get_grid_dict(location, features, settings)

            hr = _Results(super_grids=grd_dict, protein=prot)
            hr.threshold = threshold
            hr.best_island = best_island.minimal()
            hr.location = location
            hr.features = features
            hr.score = hr.score()
            hr.rank = hr._rank_features()
            return hr

    @staticmethod
    def grow_from_seed( single_grid, mask_dic, settings, prot, seed=None):
        """
        create a Extracted Hotspot object from HotspotResult object
        :param single_grid:
        :param mask_dic:
        :param settings:
        :param prot:
        :param seed:
        :return:
        """

        #single_grid.write('test.grd')
        inner = single_grid.copy_and_clear()
        inner.set_sphere(point=seed, radius=1, value=20, scaling='None')
        #mask = (sphere & single_grid) * single_grid

        #optimiser = _Results.Optimiser(mask=mask, settings=settings, peak=seed)
        #threshold, best_island = optimiser.optimize_island_threshold()
        num_gp = inner.count_grid()
        while num_gp < settings._num_gp:

            grown = Grid.grow(inner, single_grid)
            diff = grown > inner
            if diff.count_grid() < 10:
                break
            inner = grown
            num_gp = inner.count_grid()
            print(num_gp, 'out of', settings._num_gp)


        tmp_best_island = inner*single_grid
        g_vals = tmp_best_island.grid_values()
        g_vals[::-1].sort()
        try:
            threshold = g_vals[settings._num_gp]
        except IndexError:
            threshold = g_vals.min()

        best_island = grown

        print("oooh",threshold)
        if best_island is not None:
            location, features = _Results.get_interaction_type(mask_dic, best_island, threshold, settings)
            grd_dict = _Results.get_grid_dict(location, features, settings)

            hr = _Results(super_grids=grd_dict, protein=prot)
            hr.threshold = threshold
            hr.best_island = best_island.minimal()
            hr.location = location
            hr.features = features
            hr.score = hr.score()
            hr.rank = hr._rank_features()
            return hr

    @staticmethod
    def get_interaction_type(mask_dic, best_island, threshold, settings):
        """
        seperates single grid into grid by interaction type
        :return:
        """
        common_best_island = mask_dic["apolar"].common_boundaries(best_island)
        features_in_vol = {p: g * (g & common_best_island) for p, g in mask_dic.items()}
        location = features_in_vol["apolar"]
        features = _Results._get_features(interaction_dict=features_in_vol,
                                          threshold=threshold,
                                          min_feature_gp=settings.min_feature_gp)
        return location, features

    @staticmethod
    def get_grid_dict(location, features, settings):
        """
        Creates super grid dict from location and _features
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
    #     for feat in self._features:
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
    A class to handle the extraction of discrete volumes from the highest scoring regions of the hotspot maps.
    """

    class Settings(object):
        """
        Default settings for hotspot extraction

        :param float volume:
        :param float cutoff:
        :param str mode:
        :param int padding_factor:
        :param float spacing:
        :param int min_feature_gp:
        :param int max_features:
        :param float min_distance:
        :param int island_max_size:
        :param float sigma:
        :param float drug_volume:
        :param int buriedness_value:
        :param bool fragments:
        :param bool lead:
        :param bool pharmacophore:
        """

        def __init__(self, volume=150, cutoff=14, mode="seed", padding_factor=0, spacing=0.5,
                     min_feature_gp=5, max_features=10, min_distance=6, island_max_size=100, sigma=0.3,
                     drug_volume=300, buriedness_value=4, fragments=None, lead=None, pharmacophore=True):


            self.volume = volume
            self.cutoff = cutoff
            self.mode = mode
            self.padding_factor = padding_factor
            self.spacing = spacing
            self.min_feature_gp = min_feature_gp
            self.max_features = max_features
            self.min_distance = min_distance
            self.island_max_size = island_max_size
            self.sigma = sigma
            self.drug_volume = drug_volume
            self.buriedness_value = buriedness_value
            self.fragments = fragments
            self.lead = lead
            self.pharmacophore = pharmacophore

        @property
        def _num_gp(self):
            """
            number of grid point for a given volume
            :return:
            """
            return int(float(self.volume) / self.spacing ** 3)

        @property
        def _search_radius(self):
            """
            describes search radius around a given seed
            :return:
            """
            s = 3
            s += round((int(self.volume) / 50))
            print('search_radius',s)
            return s

    def __init__(self, hr, settings=None, mode="seed", volume="125", pharmacophores=True):
        """

        :param hr: An :class:`hotspots.HotspotResults` instance
        :param settings: An :class:`hotspots.Extractor.Settings` instance
        :param str mode: Options are "seed" or "global". Seed will aim to find multiple volumes, grown from suitable seed points, and works best for locating multiple hotspots in a single pocket. Global can be used to select the best target volume across the entire protein.
        :param float volume: Target volume for extractor. The extracted volume may differ slightly from the target volume
        :param bool pharmacophores: Whether to generate a pharmacophore model
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
        # if self.settings.mode == "seed":
        #     hr.super_grids = self.grid_post_process(hr.super_grids)

        # else:
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
            self._peaks = self._get_peaks()

        elif self.settings.mode == "grow":
            print("    Mode: 'grow'")
            self._peaks = self._get_peaks()

        elif self.settings.mode == "global":
            print("    Mode: 'global'")
            self._peaks = None

        else:
            raise IOError("Mode not currently supported")

        # runs and ranks extraction procedure
        print("Generate Single Grid")
        self._masked_dic, self._single_grid = Grid.get_single_grid(self.hotspot_result.super_grids)

        print("Extracting Hotspots")
        self.extracted_hotspots = self._get_extracted_hotspots()
        self._rank_extracted_hotspots()

        # generate pharmacophores
        if self.settings.pharmacophore:
            self.get_pharmacophores()

    @property
    def single_grid(self):
        return self._single_grid

    @property
    def masked_dic(self):
        return self._masked_dic

    @property
    def peaks(self):
        return self._peaks

    # def grid_post_process(self, super_grids):
    #     """
    #     carry out post-processing of fragment hotspot maps
    #
    #     Limit the size of polar islands. Keep top scores upto X grid points
    #     :return:
    #     """
    #     for probe, g in super_grids.items():
    #         if probe == "apolar":
    #             super_grids.update({probe: g.max_value_of_neighbours()})
    #
    #         else:
    #             h = g.max_value_of_neighbours()
    #             h = h.limit_island_size(self.settings.island_max_size)
    #             if h.bounding_box != super_grids["apolar"].bounding_box:
    #                 h = super_grids["apolar"].common_boundaries(g)
    #
    #             super_grids.update({probe: h})
    #
    #     # try:
    #     #     super_grids["negative"] = super_grids["negative"].deduplicate(super_grids["acceptor"],
    #     #                                                                         threshold=10,
    #     #                                                                         tolerance=2)
    #     #
    #     #     super_grids["positive"] = super_grids["positive"].deduplicate(super_grids["donor"],
    #     #                                                                         threshold=10,
    #     #                                                                         tolerance=2)
    #     # except KeyError:
    #     #     pass
    #
    #     return super_grids

    def _get_peaks(self):
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
                    peak_by_value.update({val: [(peak[0], peak[1], peak[2])]})

        average_peaks = []
        for key in peak_by_value.keys():
            x = [point[0] for point in peak_by_value[key]]
            y = [point[1] for point in peak_by_value[key]]
            z = [point[2] for point in peak_by_value[key]]
            average_peaks.append(apolar.indices_to_point(int(sum(x) / len(x)),
                                                         int(sum(y) / len(y)),
                                                         int(sum(z) / len(z))

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
            print(self.peaks)
            for peak in self.peaks:
                print(peak)
                e = _Results.from_hotspot(self.single_grid,
                                          self.masked_dic,
                                          self.settings,
                                          self.hotspot_result.protein,
                                          seed=peak)

                # if e:
                #     if e.threshold > 0:
                print(e.threshold)
                extracted_hotspots.append(e)

        elif self.settings.mode == "grow":
            print(self.peaks)
            for peak in self.peaks:
                print(peak)
                e = _Results.grow_from_seed(self.single_grid,
                                          self.masked_dic,
                                          self.settings,
                                          self.hotspot_result.protein,
                                          seed=peak)

                # if e:
                #     if e.threshold > 0:
                print(e.threshold)
                extracted_hotspots.append(e)


        else:
            e = _Results.from_hotspot(self.single_grid,
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
            -_features: islands and probes at feature point locations
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
