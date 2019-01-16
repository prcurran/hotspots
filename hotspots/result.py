#!/usr/bin/env python

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

The :mod:`hotspots.result` module contains classes to extract valuable information from the calculated Fragment Hotspot
Maps.

The main fields can be broken due using a fixed volume or fixed score mode.

This approach enables a ranking of ligand sized cavity descriptions to be
extracted. For example, these grids can be used to generate pharmacophores.

The main classes of the :mod:`hotspots.result` module are:
    -Extractor
    -Results

More information about the fragment hotspot maps method is available from:
    Radoux, C.J. et. al., Identifying the Interactions that Determine Fragment Binding at Protein Hotspots J. Med. Chem.
    2016, 59 (9), 4314-4325
    dx.doi.org/10.1021/acs.jmedchem.5b01980
"""
from __future__ import print_function, division

import copy
import operator
import tempfile
from os.path import join, dirname
from os import getcwd

import numpy as np
from ccdc.cavity import Cavity
from ccdc.io import MoleculeWriter
from ccdc.molecule import Molecule
from ccdc.protein import Protein
from scipy import optimize
from scipy.stats import percentileofscore
from skimage import feature

from hotspots.grid_extension import Grid, _GridEnsemble
from hotspots.hs_utilities import Figures
from hs_pharmacophore import PharmacophoreModel
from hs_utilities import Helper
from atomic_hotspot_calculation import AtomicHotspot, AtomicHotspotResult


class _Scorer(Helper):
    """a class to handle the annotation of objects with Fragment Hotspot Scores"""

    def __init__(self, hotspot_result, obj, tolerance):
        self.hotspot_result = hotspot_result
        self.object = obj
        self.tolerance = tolerance

        if isinstance(obj, Protein):
            self._scored_object = self.score_protein()

        elif isinstance(obj, Molecule):
            self._scored_object = self.score_molecule()

        elif isinstance(obj, Cavity):
            self._scored_object = self.score_cavity()

        elif not obj:
            self._scored_object = self.score_hotspot()

        else:
            raise IOError("supplied object not currently supported, soz!")

    @property
    def scored_object(self):
        return self._scored_object

    def score_protein(self):
        """
        score a protein's atoms, values stored as partial charge
        h_bond_distance between 1.5 - 2.5 A (2.0 A taken for simplicity)
        This method uses the cavity API to reduce the number of atoms to iterate over.
        :return:
        """
        # TODO: enable cavities to be generated from Protein objects
        feats = set([f for f in self.hotspot_result.features])
        prot = self.object
        h_bond_distance = 2.0
        interaction_pairs = {"acceptor": "donor",
                             "donor": "acceptor",
                             "pi": "apolar",
                             "aliphatic": "apolar",
                             "aromatic": "apolar",
                             "apolar": "apolar",
                             "donor_acceptor": "doneptor",
                             "dummy": "dummy"}

        cavities = Helper.cavity_from_protein(self.object)
        for cavity in cavities:
            for feature in cavity.features:
                for atm in feature.residue.atoms:
                    # score with apolar atoms
                    if atm.is_donor is False and atm.is_acceptor is False and atm.atomic_number != 1:
                        score = self.hotspot_result.super_grids['apolar'].get_near_scores(atm.coordinates)
                        if len(score) == 0:
                            score = 0
                        else:
                            score = max(score)
                        prot.atoms[atm.index].partial_charge = score

                # deal with polar atoms using cavity
                if feature.type == "acceptor" or feature.type == "donor" or feature.type =="doneptor":
                    v = feature.protein_vector
                    translate = tuple(map(h_bond_distance.__mul__, (v.x, v.y, v.z)))
                    c = feature.coordinates
                    coordinates = tuple(map(operator.add, (c.x, c.y, c.z), translate))

                    if feature.atom:
                        score = [f.score_value for f in feats if f.grid.contains_point(coordinates, tolerance=2)
                                                                and f.feature_type == interaction_pairs[feature.type]]
                        if len(score) == 0:
                            score = 0
                        else:
                            score = max(score)
                        prot.atoms[feature.atom.index].partial_charge = score

                        # score hydrogen atoms (important for GOLD)
                        a = [a.index for a in prot.atoms[feature.atom.index].neighbours
                             if int(a.atomic_number) == 1]
                        if len(a) > 0:
                            for atm in a:
                                prot.atoms[atm].partial_charge = score

        return prot

    def score_molecule(self):
        """
        score a molecule's atoms, values stored as partial charge
        :return:
        """
        # TODO: score has been placed in partial charge field. This score will persist during read and write
        mol = copy.copy(self.object)
        for atom in mol.heavy_atoms:
            atom_type = self._atom_type(atom=atom)
            score = self._score_atom_type(atom_type, atom.coordinates)
            atom.partial_charge = score

        return mol

    def score_cavity(self):
        # TODO: return scored cavity _features, the score protein function should be enough tbh
        return 0

    def score_hotspot(self, threshold=5, percentile=50):
        """
        Hotspot scored on the median value of all points included in the hotspot.
        NB: grid point with value < 5 are ommited from fragment hotspot map (hence the default threshold)
        :param percentile:
        :return:
        """
        sg = Grid.get_single_grid(self.hotspot_result.super_grids, mask=False)
        return sg.grid_score(threshold=threshold, percentile=percentile)

    # def _score_atom_type(self, grid_type, coordinates):
    #     """
    #     atom
    #     :param grid_type:
    #     :param coordinate:
    #     :param tolerance:
    #     :return:
    #     """
    #     if grid_type == "doneptor":
    #         grid_type = self._doneptor_grid(coordinates)
    #
    #     return self.hotspot_result.super_grids[grid_type].value_at_coordinate(coordinates,
    #                                                                           tolerance=self.tolerance,
    #                                                                           position=False)

    def _percentage_rank(self, obj, threshold=5):
        """
        NB: must score obj first!
        :param obj:
        :param threshold:
        :return:
        """
        mol = copy.copy(self.scored_object)
        adict = {p: g.grid_values(threshold=threshold) for p, g in self.hotspot_result.super_grids.items()}

        for atom in mol.atoms:
            atom_type = self._atom_type(atom)
            if atom_type == "doneptor":
                atom_type = self._doneptor_grid(atom.coordinates)
            atom.partial_charge = percentileofscore(adict[atom_type], atom.partial_charge)

        return mol

    # def _doneptor_grid(self, coordinates):
    #     """
    #     An atom is scored from the grid which yields the highest value
    #     :param coordinates:
    #     :param grid_type:
    #     :return:
    #     """
    #     scores = [self.hotspot_result.super_grids["donor"].value_at_coordinate(coordinates,
    #                                                                            tolerance=self.tolerance,
    #                                                                            position=False),
    #               self.hotspot_result.super_grids["acceptor"].value_at_coordinate(coordinates,
    #                                                                               tolerance=self.tolerance,
    #                                                                               position=False)
    #               ]
    #     d = dict(zip(scores, ["donor", "acceptor"]))
    #     return d[max(d.keys())]

    @staticmethod
    def _atom_type(atom):
        """
        from a ccdc Atom, the "atom type" is returned
        :param a:
        :return:
        """
        if atom.is_donor and atom.is_acceptor:
            return "doneptor"

        elif atom.is_acceptor:
            return "acceptor"

        elif atom.is_donor:
            return "donor"

        elif atom.atomic_symbol == "Xe":
            return "dummy"

        else:
            return "apolar"


class Results(object):
    """
    A Hotspot_results object is returned at the end of a Hotspots calculation. It contains functions for accessing
    and using the results.
    """

    def __init__(self, super_grids, protein, buriedness=None, pharmacophore=None):
        try:
            self.super_grids = super_grids
            for probe, g in super_grids.items():
                b = g.bounding_box
        except:
            raise TypeError("Not a valid Grid")

        self.protein = protein
        self.buriedness = buriedness
        self.pharmacophore = None
        self.features = self._get_features(interaction_dict=super_grids)

        if pharmacophore:
            self.pharmacophore = self.get_pharmacophore_model()

    class _HotspotFeature(object):
        """
        class to hold polar islands above threshold "_features"
        purpose: enables feature ranking
        """

        def __init__(self, feature_type, grid):
            """

            :param feature_type:
            :param grid:
            """
            self._feature_type = feature_type
            self._grid = grid
            self._feature_coordinates = grid.centroid()
            self._count = (grid > 0).count_grid()
            self._score_value = self.score_feature()
            self._sphere = None

            # set these
            self._rank = None
            self.superstar_results = []

        @property
        def feature_type(self):
            return self._feature_type

        @property
        def grid(self):
            return self._grid

        @property
        def feature_coordinates(self):
            return self._feature_coordinates

        @property
        def sphere(self):
            return self._sphere

        @property
        def count(self):
            return self._count

        @property
        def score_value(self):
            return self._score_value

        @property
        def rank(self):
            return self._rank

        def score_feature(self):
            """
            returns
            :return:
            """
            nx, ny, nz = self.grid.nsteps

            return sum([self.grid.value(i, j, k)
                        for i in range(nx) for j in range(ny) for k in range(nz)
                        if self.grid.value(i, j, k) > 0]) / self.count

    def tractability_scores(self):
        extractor_settings = Extractor.Settings()
        extractor_settings.cutoff = 5
        extractor_settings.island_max_size = 500

        extractor = Extractor(self, mode="global", volume=500, settings=extractor_settings)
        hist = extractor.extracted_hotspots[0].get_map_values()
        all = []
        for x in hist.values():
            all += x.tolist()

        return np.median(all)


    def score(self, obj=None, tolerance=2):
        """
        Given a supported CCDC object, will return the object annotated with Fragment Hotspot scores

        :param obj:
        :param return_value:
        :return:

        TODO: Complete this docstring
        """

        return _Scorer(self, obj, tolerance).scored_object

    def _filter_map(self, g1, g2, tol):
        """
        Experimental feature.
        Takes 2 grids of the same size and coordinate frames. Points that are
        zero in one grid but sampled in the other are
        set to the mean of their nonzero neighbours. Grids are then subtracted and
        the result is returned.
        :param int tol: how many grid points away to consider scores from
        :param g1: a :class: "ccdc.utilities.Grid" instance
        :param g2: a :class: "ccdc.utilities.Grid" instance
        :return: a :class: "ccdc.utilities.Grid" instance
        """

        def filter_point(x, y, z):
            loc_arr = np.array(
                [g[x + i][y + j][z + k] for i in range(-tol, tol + 1) for j in range(-tol, tol + 1) for k in
                 range(-tol, tol + 1)])
            if loc_arr[loc_arr > 0].size != 0:
                #print(np.mean(loc_arr[loc_arr > 0]))
                new_grid[x][y][z] = np.mean(loc_arr[loc_arr > 0])

        vfilter_point = np.vectorize(filter_point)
        com_bound_box = g1.bounding_box
        com_spacing = g1.spacing

        arr1 = g1.get_array()
        arr2 = g2.get_array()

        b_arr1 = np.copy(arr1)
        b_arr2 = np.copy(arr2)

        b_arr1[b_arr1 > 0] = 1.0
        b_arr2[b_arr2 > 0] = -1.0

        diff_arr = b_arr1 + b_arr2

        unmatch1 = np.where(diff_arr == 1)
        unmatch2 = np.where(diff_arr == -1)

        g = arr1
        new_grid = np.copy(arr1)
        vfilter_point(unmatch2[0], unmatch2[1], unmatch2[2])
        f_arr1 = np.copy(new_grid)

        g = arr2
        new_grid = np.copy(arr2)
        vfilter_point(unmatch1[0], unmatch1[1], unmatch1[2])
        f_arr2 = np.copy(new_grid)

        sel_arr = f_arr1 - f_arr2
        sel_arr[sel_arr < 0] = 0
        sel_map = Grid(origin=com_bound_box[0], far_corner=com_bound_box[1], spacing=com_spacing, _grid=None)

        idxs = sel_arr.nonzero()
        vals = sel_arr[idxs]

        as_triads = zip(*idxs)
        for (i, j, k), v in zip(as_triads, vals):
            sel_map._grid.set_value(int(i), int(j), int(k), v)

        return sel_map

    def get_difference_map(self, other, tolerance):
        """
        Experimental feature.
        Generates maps to highlight selectivity for a target over an off target cavity. Proteins should be aligned
        by the binding site of interest prior to calculation.
        High scoring regions of a map represent areas of favourable interaction in the target binding site, not
        present in off target binding site
        :param other: a :class:`hotspots.result.Results` instance
        :param int tolerance: how many grid points away to apply filter to
        :return: a :class:`hotspots.result.Results` instance
        """

        selectivity_grids = {}
        for probe in self.super_grids.keys():
            g1 = self.super_grids[probe]
            g2 = other.super_grids[probe]
            og1, og2 = Grid.common_grid([g1, g2])
            sele = self._filter_map(og1, og2, tolerance)
            selectivity_grids[probe] = sele
        hr = Results(selectivity_grids, self.protein, None, None)
        return hr

    @staticmethod
    def from_grid_ensembles(res_list, prot_name, charged=False):
        """
        Experimental feature.
        Creates ensemble map from a list of Results. Structures in the ensemble have to aligned by the
        binding site of interest prior to the hotspots calculation.
        :param res_list: list of hotspots.result.Results
        :param str prot_name: str
        :param str out_dir: str
        :return: a :class:`hotspots.result.Results` instance
        """
        if charged:
            probe_list = ["acceptor", "apolar", "donor", "positive", "negative"]
        else:
            probe_list = ["acceptor", "apolar", "donor"]

        grid_dic = {}

        for p in probe_list:
            grid_list_p = [r.super_grids[p].minimal() for r in res_list]
            ens = _GridEnsemble()
            grid_dic[p] = ens.from_grid_list(grid_list_p, getcwd(), prot_name, p)

        hr = Results(grid_dic, protein=res_list[0].protein)
        return hr

    def get_pharmacophore_model(self, identifier="id_01", cutoff=5):
        """
        Generates a :class:`hotspots.hotspot_pharmacophore.PharmacophoreModel` instance from peaks in the hotspot maps


        :param str identifier: Identifier for displaying multiple models at once
        :param float cutoff: The score cutoff used to identify islands in the maps. One peak will be identified per island
        :return: a :class:`hotspots.hotspot_pharmacophore.PharmacophoreModel` instance
        """
        return PharmacophoreModel.from_hotspot(self.protein, self.super_grids, identifier=identifier, cutoff=cutoff)

    def get_map_values(self):
        """
        get the number zero grid points for the Fragment Hotspot Result
        :return: dict of str(probe type) by a :class:`numpy.array` (non-zero grid point scores)
        """
        data = Figures.histogram(self, False)
        return data

    def get_histogram(self, fpath="histogram.png"):
        """
        get histogram of zero grid points for the Fragment Hotspot Result
        :param fpath: path to output file
        :return:
        """
        data, plt = Figures.histogram(self, plot=True)
        plt.savefig(fpath)
        return data, plt

    # def get_2D_diagram(self, ligand, fpath="diagram.png", title=False):
    #     """
    #
    #     :param ligand:
    #     :param fpath:
    #     :param title:
    #     :return:
    #     """
    #     Figures._2D_diagram(hr, ligand, title=False, output="diagram.png")

    def get_elaboration_potential(self, large_cavities):
        """
        Is the hotspot within a drug size cavity:
        0 = False
        1 = True
        TODO: develop a more sophicated method to evaluate elaboration potential
        TODO: Complete docstring
        :param large_cavities:
        :return:
        """
        centroid = self.location.centroid()
        cavity = [c for c in large_cavities if c.contains_point(centroid)]

        if len(cavity) == 0:
            return 0
        else:
            return 1

    def _get_superstar_profile(self, feature_radius=1.5, nthreads=6, features=None, best_volume=None):
        """
        EXPERIMENTAL
        :return:
        """
        # set additional object properties
        if features:
            self.features = features
        else:
            self.features = self._get_features(threshold=5, min_feature_size=6)

        if best_volume:
            self.best_volume = best_volume
        else:
            self.best_volume = Grid.get_single_grid(self.super_grids, mask=False)

        self.feature_spheres = self.best_volume.copy_and_clear()
        for feat in self.features:
            self.feature_spheres.set_sphere(point=feat.feature_coordinates,
                                            radius=feature_radius,
                                            value=1,
                                            scaling="None"
                                            )

        # superstar run
        centroid = [self.best_volume.centroid()]
        a = AtomicHotspot()
        a.settings.atomic_probes = ["carbonyl_oxygen", "carboxylate", "pyramidal_r3n", "water_oxygen"]

        self.superstar_result = a.calculate(protein=self.protein,
                                            nthreads=nthreads,
                                            cavity_origins=centroid)

        self.ss = []

        # find overlap
        for r in self.superstar_result:
            common_spheres, common_result = Grid.common_grid([self.feature_spheres, r.grid])
            r.grid = (common_spheres & common_result) * common_result

        # assign island to Hotspot Feature
        feat_id = []
        ss_id = []
        score = []
        import pandas as pd

        for i, feat in enumerate(self.features):

            for r in self.superstar_result:
                feat_id.append(i)
                ss_id.append(r.identifier)

                ss_dict = {Helper.get_distance(feat.feature_coordinates, island.centroid()): island
                           for island in r.grid.islands(threshold=1)
                           if Helper.get_distance(feat.feature_coordinates, island.centroid()) < 1}

                if len(ss_dict) == 0:
                    g = r.grid.copy_and_clear()

                else:
                    shortest = sorted([f[0] for f in ss_dict.items()], reverse=False)[0]
                    g = ss_dict[shortest]

                feat.superstar_results.append(AtomicHotspotResult(identifier=r.identifier,
                                                                  grid=g,
                                                                  buriedness=None)
                                              )

                score.append(g.grid_score(threshold=1, percentile=50))

        return pd.DataFrame({"feature_id": feat_id, "interaction": ss_id, "score": score})

    @staticmethod
    def _get_features(interaction_dict, threshold=5, min_feature_gp=6, excluded=("apolar")):
        """
        returns Hotspot Feature object with a score to enable ranking
        :param probe:
        :param g:
        :return:
        """
        f = []
        for probe, g in interaction_dict.items():
            if len(g.islands(threshold=threshold)) > 0:
                for island in g.islands(threshold=threshold):
                    if (island > threshold).count_grid() > min_feature_gp and probe not in excluded:
                        f.append(Results._HotspotFeature(probe, island))
        return f

    def _rank_features(self):
        """
        rank _features based upon feature score (TO DO: modify score if required)
        :return:
        """
        feature_by_score = {feat.score_value: feat for feat in self.features}
        score = sorted([f[0] for f in feature_by_score.items()], reverse=True)
        for i, key in enumerate(score):
            feature_by_score[key]._rank = int(i + 1)

    def extract_pocket(self, whole_residues=False):
        """
        Create a :class:`ccdc.Protein` containing atoms or residues that have a score > 0

        :param bool whole_residues: whether to include all residue atoms if only a subset have a score > 0
        :return: a :class:`ccdc.Protein` instance
        """
        prot_scores = self.score_protein()
        pocket = self.protein.copy()
        pocket.remove_hydrogens()
        for residue in pocket.residues:
            keep_residue = False
            for atom in residue.atoms:
                a_id = "{0}/{1}/{2}".format(residue.chain_identifier, residue.identifier.split(':')[1][3:],
                                            atom.label)
                atom_type = self._get_atom_type(atom)

                if atom_type == 'doneptor':
                    score = max([prot_scores[a_id]['donor'], prot_scores[a_id]['acceptor']])
                else:
                    score = prot_scores[a_id][atom_type]

                if score > 0:
                    keep_residue = True
                elif score == 0 and not whole_residues:
                    pocket.remove_atom(atom)

            if whole_residues and not keep_residue:
                pocket.remove_atoms(residue.atoms)
        return pocket

    def _ngl_widget(self, out_dir=None):
        """
        jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000000.0
        creates ngl widget from hotspot. For use in ipython notebooks
        :param str out_dir:
        :return:
        """
        import nglview as nv
        from ipywidgets import IntSlider, interact

        color_dict = {"apolar": "yellow",
                      "donor": "blue",
                      "acceptor": "red",
                      "negative": "magenta",
                      "positive": "cyan"}
        if out_dir:
            out = Helper.get_out_dir(out_dir)
        else:
            out = tempfile.mkdtemp()

        for p, g in self.super_grids.items():
            g.write(join(out, "{}.ccp4".format(p)))

        with MoleculeWriter(join(out, "protein.pdb")) as w:
            w.write(self.protein)

        view = nv.NGLWidget()
        view.add_component(join(out, "protein.pdb"))

        k = self.super_grids.keys()
        for i, p in enumerate(k):
            view.add_component(join(out, "{}.ccp4".format(p)))
            view.add_representation('isosurface', component=i + 1)
            view.update_representation(component=i + 1, color=color_dict[p])

        @interact(x=IntSlider(description="HS Score", min=0, max=30, step=1))
        def f(x):
            view.update_representation(component=1, isolevel=int(x), isoleveltype='value')
            view.update_representation(component=2, isolevel=int(x), isoleveltype='value')
            view.update_representation(component=3, isolevel=int(x), isoleveltype='value')
            view.update_representation(component=4, isolevel=int(x), isoleveltype='value')
            view.update_representation(component=5, isolevel=int(x), isoleveltype='value')

        return view


class _ExtractedHotspot(object):
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
            # threshold = optimize.fminbound(self._count_island_points, 0, 50, xtol=0.025)
            x0 = np.array([14])
            b = optimize.Bounds(0, 50)
            ret = optimize.minimize_scalar(self._count_island_points, x0, bounds=(0, 50), method='bounded')
            threshold = ret.x
            if threshold > 48:
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

    def __init__(self, hr):
        self.result = hr

    def from_hotspot(self, single_grid, mask_dic, settings, prot, seed=None):
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

        optimiser = _ExtractedHotspot.Optimiser(mask=mask, settings=settings, peak=seed)
        threshold, best_island = optimiser.optimize_island_threshold()

        if best_island is not None:
            location, features = self.get_interaction_type(mask_dic, best_island, threshold, settings)
            grd_dict = self.get_grid_dict(location, features, settings)

            hr = Results(super_grids=grd_dict, protein=prot)


            hr.threshold = threshold
            hr.best_island = best_island.minimal()
            hr.location = location
            hr.features = features
            hr.score_value = hr.score()
            hr.rank = hr._rank_features()
            return hr

    def grow_from_seed(self, single_grid, mask_dic, settings, prot, seed=None):
        """
        create a Extracted Hotspot object from HotspotResult object
        :param single_grid:
        :param mask_dic:
        :param settings:
        :param prot:
        :param seed:
        :return:
        """

        inner = single_grid.copy_and_clear()
        inner.set_sphere(point=seed, radius=1, value=20, scaling='None')
        # mask = (sphere & single_grid) * single_grid

        # optimiser = _Results.Optimiser(mask=mask, settings=settings, peak=seed)
        # threshold, best_island = optimiser.optimize_island_threshold()
        num_gp = inner.count_grid()
        grown = Grid.grow(inner, single_grid)
        while num_gp < settings._num_gp:

            grown = Grid.grow(inner, single_grid)
            diff = grown > inner
            if diff.count_grid() < 10:
                break
            inner = grown
            num_gp = inner.count_grid()
            print(num_gp, 'out of', settings._num_gp)

        tmp_best_island = inner * single_grid
        g_vals = tmp_best_island.grid_values()
        g_vals[::-1].sort()
        try:
            threshold = g_vals[settings._num_gp]
        except IndexError:
            threshold = g_vals.min()

        best_island = grown

        if best_island is not None:
            location, features = _ExtractedHotspot.get_interaction_type(mask_dic, best_island, threshold, settings)
            grd_dict = _ExtractedHotspot.get_grid_dict(location, features, settings)

            hr = Results(super_grids=grd_dict, protein=prot)

            hr.threshold = threshold
            hr.best_island = best_island.minimal()
            hr.location = location
            hr.features = features
            hr.score_value = hr.score()
            hr.rank = hr._rank_features()
            return hr

    def get_interaction_type(self, mask_dic, best_island, threshold, settings):
        """
        seperates single grid into grid by interaction type
        :return:
        """
        common_best_island = mask_dic["apolar"].common_boundaries(best_island)
        features_in_vol = {p: g * (g & common_best_island) for p, g in mask_dic.items()}
        location = features_in_vol["apolar"]
        features = self.result._get_features(interaction_dict=features_in_vol,
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
        feature_by_score = {f.score_value: f for f in features}
        features = [feature_by_score[s]
                    for s in sorted([f[0] for f in feature_by_score.items()], reverse=True)][:settings.max_features - 1]
        for probe in interaction_types:
            if settings.mode == "seed":
                grids = [feat.grid for feat in features
                         if feat.feature_type == probe and
                         feat.score_value >= settings.cutoff]

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
            print('search_radius', s)
            return s

    def __init__(self, hr, settings=None, mode="seed", volume="125", pharmacophores=True, mvon=True):
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
        self._single_grid = None
        self._masked_dic = None
        self.settings.mode = mode
        self.settings.volume = volume
        self.settings.pharmacophore = pharmacophores
        self.out_dir = None
        self.mvon = mvon

        # fragment hotspot post processing
        # if self.settings.mode == "seed":
        #     hr.super_grids = self.grid_post_process(hr.super_grids)

        # else:
        if self.mvon:
            hr.super_grids.update({probe: g.max_value_of_neighbours() for probe, g in hr.super_grids.items()})

        try:
            hr.super_grids["negative"] = hr.super_grids["negative"].deduplicate(hr.super_grids["acceptor"],
                                                                                threshold=10,
                                                                                tolerance=2)

            hr.super_grids["positive"] = hr.super_grids["positive"].deduplicate(hr.super_grids["donor"],
                                                                                threshold=10,
                                                                                tolerance=2)
        except KeyError:
            print("deduplication failed, make sure grids have common boundaries")
            pass

        hr.super_grids.update({probe: g.minimal() for probe, g in hr.super_grids.items()})

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

                extracted = _ExtractedHotspot(self.hotspot_result)
                e = extracted.from_hotspot(self.single_grid,
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
                extracted = _ExtractedHotspot(self.hotspot_result)
                e = extracted.grow_from_seed(self.single_grid,
                                             self.masked_dic,
                                             self.settings,
                                             self.hotspot_result.protein,
                                             seed=peak)

                # if e:
                #     if e.threshold > 0:
                print(e.threshold)
                extracted_hotspots.append(e)


        else:
            extracted = _ExtractedHotspot(self.hotspot_result)
            e = extracted.from_hotspot(self.single_grid,
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
        hotspot_by_score = {hotspot.score_value: hotspot for hotspot in self.extracted_hotspots}
        score = sorted([f[0] for f in hotspot_by_score.items()], reverse=True)

        for i, key in enumerate(score):
            hotspot_by_score[key].rank = int(i + 1)

        extracted_hotspots_by_rank = {h.rank: h for h in self.extracted_hotspots}
        self.extracted_hotspots = [value for (key, value) in sorted(extracted_hotspots_by_rank.items())]

        for i, hs in enumerate(self.extracted_hotspots):
            hs.identifier = "rank_{}".format(hs.rank)
            print("rank", hs.rank, "score", hs.score_value)

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
