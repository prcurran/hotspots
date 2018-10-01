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

The main classes of the :mod:`fragment_hotspot_maps.extraction` module are:
    -BuildLocation
    -BuildFeature
    -HotspotBuilder
"""
from __future__ import print_function, division
from os.path import join, dirname, exists
import numpy as np
import numpy.ma as ma

from ccdc.io import MoleculeWriter
from utilities import Utilities
from enhanced_grid import Grid

from skimage import feature
from os import mkdir


class BuildLocation(Utilities):
    """
    class to hold information about apolar interactions
    """
    def __init__(self, apolar, indices, identifier, **kwargs):
        self.identifier = identifier
        self._apolar = apolar
        self._indices = indices
        self._kwargs = kwargs

        self.assigned_features = []
        self._peak_score = apolar.value(int(indices[0]), int(indices[1]), int(indices[2]))
        self._coordinates = apolar.indices_to_point(indices[0], indices[1], indices[2])
        self._parent_island = self.get_parent_island()
        self._island = self.get_surrounding_points()

        # set by HotspotBuilder class
        self.extracted_super_grids = None
        self.extracted_probes = None

    @property
    def apolar(self):
        """apolar grid from fragment hotspot maps calculation (gaussian filter applied)"""
        return self._apolar

    @property
    def indices(self):
        """indices of local maxima within the apolar grid"""
        return self._indices

    @property
    def peak_score(self):
        """grid score of local maxima within the apolar grid"""
        return self._peak_score

    @property
    def coordinates(self):
        """coordinates of local maxima within the apolar grid"""
        return self._coordinates

    @property
    def settings(self):
        """settings required for `fragment_hotspot_maps.extraction.BuildLocation` class"""
        _settings = HotspotBuilder.Settings(self._kwargs)
        _settings.mode = self._kwargs.get("mode")
        return _settings

    @property
    def parent_island(self):
        """island (@settings.cutoff) of the apolar grid in which the local maxima is located"""
        return self._parent_island

    @property
    def island(self):
        """extracted island from either volume or score mode"""
        return self._island

    def get_parent_island(self):
        """
        find apolar island(@settings.cutoff) which contains local maxima
        :return: self._parent_island
        """
        if self.settings.mode == "volume":
            islands = self._apolar.islands(threshold=self.settings.cutoff - 2)

        elif self.settings.mode == "score":
            islands = self._apolar.islands(threshold=17)

        else:
            raise RuntimeError("Not a valid extraction mode")

        for island in islands:
            if island.contains_point(self.coordinates, tolerance=2):
                self._parent_island = island
                return self._parent_island
            else:
                continue

        raise RuntimeError("No parent island found with valid peak coordinates")

    def weight_by_distance(self):
        """
        modification of parent island, grid points are weighted by distance from local maxima
        :return: adjusted_parent
        """
        adjusted_parent = self.parent_island.copy_and_clear()
        nx, ny, nz = self.parent_island.nsteps
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if self.parent_island.value(i, j, k) < (self.settings.cutoff - 5):
                        adjusted_parent.set_value(i, j, k, 0)
                    else:
                        coords = self.parent_island.indices_to_point(i, j, k)
                        dist = self._get_distance(coords, self.coordinates)
                        new_value = (1 / (0.1 + dist)) * self.parent_island.value(i, j, k)
                        adjusted_parent.set_value(i, j, k, new_value)
        return adjusted_parent

    def get_surrounding_points(self):
        """
        points surround the local maxima are collected:
            -volume mode = a user defined volume of grid points is selected from weighted parent island
            (default = 75A^3).
            -score mode = a user defined score is used to select grid points
            (default = 17).
        :return: self._island
        """
        if self.settings.mode == "volume":
            print("Fixed volume mode")
            adjusted_parent = self.weight_by_distance()
            fragment_volume = adjusted_parent.restricted_volume(self.settings.volume)
            self._island = fragment_volume.gaussian(0.3)
            return self._island

        elif self.settings.mode == "score":
            print("Score cutoff mode")
            maxima = 400 / 0.125            # max size 400 A^3
            array = ma.masked_less_equal(self.parent_island.get_array(), 17.0).compressed()
            if len(array) > maxima:
                vol = self.parent_island.restricted_volume(volume=400)
            else:
                vol = self.parent_island
            self._island = vol.gaussian(0.2)
            return self._island

        else:
            raise RuntimeError("Not a valid extraction mode")


class BuildFeature(object):
    """
    class to hold information about polar interactions points
    """
    def __init__(self, identifier, feature_type, island, grid, probes):
        self.identifier = identifier
        self._feature_type = feature_type
        self._island = Grid.super_grid(2, island)
        self._grid = grid
        self._probes = probes

        # generate data
        self._coordinates = island.centroid()
        self._island_probes = self.get_island_probes()
        self._island_maxima = island.extrema[1]

    @property
    def feature_type(self):
        """
        feature type is derived from the input grid.
        currently supports:
            -donor
            -acceptor
            -negative
            -positive
        """
        return self._feature_type

    @property
    def island(self):
        """feature island (@settings.cutoff)"""
        return self._island

    @property
    def grid(self):
        """polar grid from fragment hotspot maps calculation"""
        return self._grid

    @property
    def probes(self):
        """all probes for feature interaction type"""
        return self._probes

    @property
    def coordinates(self):
        """coordinates of feature centroid"""
        return self._coordinates

    @property
    def island_probes(self):
        """sampled probes which have their priority atom within the feature island bounding box"""
        return self._island_probes

    @property
    def island_maxima(self):
        """value of the maximum score in the island"""
        return self._island_maxima

    def get_island_probes(self):
        """
        assign probes to island, return self._island_probes
        :return: self._island_probes
        """
        island_probes = []
        for mol in self._probes:
            coords = mol.atoms[0].coordinates
            if self._island.contains_point(coords, tolerance=2):
                island_probes.append(mol)
        self._island_probes = island_probes
        return self._island_probes


class HotspotBuilder(Utilities):
    """
    A class to handle the extraction of discrete, fragment size hotspots from the original maps
    """

    def __init__(self, super_grids, sampled_probes, out_dir, protein, kw):
        self.super_grids = {}
        for probe, g in super_grids.items():
            print(probe)
            if probe == "apolar":
                self.super_grids.update({probe: g.gaussian(0.2)})
            else:
                self.super_grids.update({probe: g.gaussian(0.5)})

        self.super_grids["negative"] = self.super_grids["negative"].deduplicate(self.super_grids["acceptor"],
                                                                                threshold=14,
                                                                                tolerance=0)
        self.super_grids["positive"] = self.super_grids["positive"].deduplicate(self.super_grids["donor"],
                                                                                threshold=14,
                                                                                tolerance=0)
        self.sampled_probes = sampled_probes
        self.out_dir = out_dir
        self.protein = protein
        self.settings = self.Settings(kw)

        self._locations = self.get_locations()
        self._features = self.get_features()
        self.format_data()

    class Settings(object):
        """
        Default settings for hotspot extraction
        """

        def __init__(self, kw):
            self.cutoff = kw.get("cutoff", 14)
            self.volume = kw.get("volume", 75)
            self.grid_points = int(float(self.volume) / 0.125)
            self.min_grid_points = 100
            self.max_probes = 50
            self.mode = None
            self.distance_cutoff = 8
            self.match_threshold = 0.2

    @property
    def extracted_super_grids(self):
        """
        grid dictionary formatted to create a:
        `fragment_hotspot_maps.fragment_hotspot_maps.HotspotResult` class instance
        """
        return self._extracted_super_grids

    @property
    def extracted_probes(self):
        """
        probe dictionary formatted to create a:
        `fragment_hotspot_maps.fragment_hotspot_maps.HotspotResult` class instance
        """
        return self._extracted_probes

    @property
    def locations(self):
        """list of `fragment_hotspot_maps.extraction.BuildLocation` class instances"""
        return self._locations

    @property
    def features(self):
        """list of unassigned `fragment_hotspot_maps.extraction.BuildFeature` class instances"""
        return self._features

    def get_locations(self):
        """
        locate peaks in apolar maps and define fragment size volume
        :return: self._locations
        """
        locations = []
        apolar = self.super_grids["apolar"]
        apolar_array = apolar.get_array()
        peaks = feature.peak_local_max(apolar_array, min_distance=6)
        for i, peak in enumerate(peaks):
            if apolar.value(int(peak[0]), int(peak[1]), int(peak[2])) > self.settings.cutoff:
                try:
                    build_locations = BuildLocation(apolar, peak, i, mode="volume")
                    locations.append(build_locations)
                except RuntimeError:
                    pass
            else:
                continue
        self._locations = locations
        return self._locations

    def get_features(self):
        """
        generate InteractionFeatures which contains all the information required for hotspot feature assignment
        :return:
        """
        features = []
        for probe, grid in self.super_grids.items():
            for i, island in enumerate(grid.islands(threshold=12)):
                # filter noise
                points = len(ma.masked_less_equal(island.get_array(), self.settings.cutoff).compressed())
                if points > 10 and probe != "apolar":
                    features.append(BuildFeature(identifier="{}_{}".format(probe, i),
                                                 feature_type=probe,
                                                 island=island,
                                                 grid=grid,
                                                 probes=self.sampled_probes[probe]
                                                 )
                                    )
                else:
                    continue
        self._features = features
        return self._features

    def construct(self):
        """
        assigns hotspot features to locations. attached to `fragment_hotspot_maps.extraction.BuildLocation`
        class instance.
        """
        for feat in self.features:
            # 1) distance filter
            distances = {self._get_distance(location.coordinates, feat.coordinates): location
                         for location in self.locations
                         if self._get_distance(location.coordinates, feat.coordinates) <
                         self.settings.distance_cutoff}

            # 2) probe location filter
            num_locations = len(distances.values())
            match_score = {}
            if num_locations != 0:
                for i, location in enumerate(distances.values()):
                    all_probes = len(feat.island_probes)
                    assigned_probes = [mol for mol in feat.island_probes
                                       if location.island.contains_point(mol.atoms[3].coordinates)]
                    match_score.update({(len(assigned_probes) / all_probes): location})
                location_key = sorted(match_score.items(), key=lambda x: x[0], reverse=True)[0][0]
                location = match_score[location_key]
                location.assigned_features.append(feat)
            else:
                print("Unable to assign feature")

    def format_data(self):
        """
        create grid_dic and probe_dic to create hotspot results object
        """
        for location in self.locations:
            grid_dic = {}
            probe_dic = {}
            for probe in self.super_grids.keys():
                if probe == "apolar":
                    grid_dic.update({probe: location.island})
                else:
                    interaction_islands = [feat.island for feat in location.assigned_features
                                           if probe == feat.feature_type]
                    interaction_probes = [mol for feat in location.assigned_features
                                          for mol in feat.island_probes
                                          if probe == feat.feature_type]
                    if len(interaction_islands) > 0:
                        grid = Grid.super_grid(1, *interaction_islands)
                    else:
                        grid = location.island.copy_and_clear()
                    grid_dic.update({probe: grid})
                    probe_dic.update({probe: interaction_probes})
                location.extracted_super_grids = grid_dic
                location.extracted_probes = probe_dic

    def write(self, **kwargs):
        """
        write out information to aid debugging: valid modes:
            -locations: spheres and islands at apolar peak locations
            -features: islands and probes at feature point locations
        """
        mode = kwargs.get("mode")
        pymol_out = ""
        pymol_out += 'from pymol import cmd\nfrom pymol.cgo import *\n'

        if mode == "locations":
            for i, location in enumerate(self.locations):
                od = join(self.out_dir, "apolar_islands")
                if not exists(od):
                    mkdir(od)
                location.island.write(join(od, "id_{}.grd".format(i)))

                sphere = 'id_{0} = [COLOR, 1.00, 1.000, 0.000] + ' \
                         '[ALPHA, 0.8] + ' \
                         '[SPHERE, float({1}), float({2}), float({3}), float(0.5)]\n' \
                    .format(location.identifier,
                            location.coordinates[0],
                            location.coordinates[1],
                            location.coordinates[2]
                            )

                pymol_out += sphere
                pymol_out += '\ncmd.load_cgo(id_{0}, "id_{0}", 1)' \
                    .format(location.identifier)
                pymol_out += '\ncmd.group("Peaks", members="id_{0}")\n' \
                    .format(location.identifier)

                pymol_out += """
nh = {0}
for n in range(nh):
    cmd.load(r'apolar_islands/id_%s.grd' % (n), 'apolar_%s' % (n))
    cmd.isosurface('surface_apolar_%s' % (n), 'apolar_%s' % (n), 5)
    cmd.set('transparency', 0.7, 'surface_apolar_%s' % (n))
    cmd.color('yellow', 'surface_apolar_%s' % (n))
for n in range(nh):
    cmd.group('hotspot_%s'%(n), members= 'surface_apolar_%s'%(n))
    cmd.group('hotspot_%s'%(n), members= 'apolar_%s'%(n))""" \
                    .format(len(self.locations))

        elif mode == "features":
            pymol_out = ""
            pymol_out += 'from pymol import cmd\nfrom pymol.cgo import *\n'
            pymol_out += """
colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'br4', 'positive':'cyan'}"""
            for j, feat in enumerate(self.features):
                od = join(self.out_dir, "features")
                if not exists(od):
                    mkdir(od)

                island_fname = "{}_{}.grd".format(feat.feature_type, j)
                probe_fname = "{}_{}_probes.mol2".format(feat.feature_type, j)

                feat.island.write(join(od, island_fname))

                with MoleculeWriter(join(od, probe_fname)) as writer:
                    for mol in feat.island_probes:
                        writer.write(mol)

                pymol_out += """
#cmd.load("features/{0}",'probe_{1}')
cmd.load("features/{2}",'{3}_{1}')
cmd.isosurface('surface_{3}_{1}', '{3}_{1}', 5)
cmd.set('transparency', 0.7, 'surface_{3}_{1}')
cmd.color(colour_dict['{3}'], 'surface_{3}_{1}')
#cmd.group('feature_{1}', members= 'surface_{3}_{1}')
#cmd.group('feature_{1}', members= 'probe_{1}')
                """. format(probe_fname, j, island_fname, feat.feature_type)

        else:
            raise RuntimeError("write mode not supported")

        with open(join(self.out_dir, "{}.py".format(mode)), "w") as pymol_file:
            pymol_file.write(pymol_out)
