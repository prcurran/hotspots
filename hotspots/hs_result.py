#!/usr/bin/env python

"""
More information about the fragment hotspot maps method is available from:
    Radoux, C.J. et. al., Identifying the Interactions that Determine Fragment Binding at Protein Hotspots J. Med. Chem.
    2016, 59 (9), 4314-4325
    dx.doi.org/10.1021/acs.jmedchem.5b01980
"""
from __future__ import print_function, division

import copy
import multiprocessing
import numpy as np
import operator
import pkg_resources
import random
import sys
import tempfile
import time
from ccdc.cavity import Cavity
from ccdc.io import MoleculeWriter, MoleculeReader
from ccdc.molecule import Molecule, Coordinates
from ccdc.protein import Protein
from ccdc.utilities import PushDir
from os import system, environ
from os.path import join
# from best_volume import Extractor
from scipy.stats import percentileofscore
from tqdm import tqdm

#from atomic_hotspot_calculation import AtomicHotspot, AtomicHotspotResult
from grid_extension import Grid
#from hotspots.atomic_hotspot_calculation import AtomicHotspot, AtomicHotspotResult
from hotspots.grid_extension import Grid
from hotspots.hs_pharmacophore import PharmacophoreModel
from hotspots.hs_utilities import Figures, Helper
from hs_pharmacophore import PharmacophoreModel


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

    def extract_volume(self):
        pass

    def score(self, obj=None, tolerance=2):
        """
        Given a supported CCDC object, will return the object annotated with Fragment Hotspot scores

        :param obj:
        :param return_value:
        :return:

        TODO: Complete this docstring
        """

        return _Scorer(self, obj, tolerance).scored_object

    def get_selectivity_map(self, other):
        """
        Generate maps to highlight selectivity for a target over an off target cavity. Proteins should be aligned
        by the binding site of interest prior to calculation.

        High scoring regions of a map represent areas of favourable interaction in the target binding site, not
        present in off target binding site

        :param other: a :class:`hotspots.hotspot_calculation.HotspotResults` instance
        :return: a :class:`hotspots.hotspot_calculation.HotspotResults` instance
        """

        selectivity_grids = {}
        for probe in self.super_grids.keys():
            g1 = self.super_grids[probe]
            g2 = other.super_grids[probe]
            og1, og2 = self._common_grid(g1, g2)
            sele = og1 - og2
            selectivity_grids[probe] = sele
        hr = Runner.HotspotResults(selectivity_grids, self.protein, self.fname, None, None, self.out_dir)
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
        data, plt = Figures.histogram(self, plot)
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
    def _get_features(interaction_dict, threshold=5, min_feature_gp=6, excluded=["apolar"]):
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


class _Scorer(object):
    """a class to handle the annotation of objects with Fragment Hotspot Scores"""

    def __init__(self, hotspot_result, object, tolerance):
        self.hotspot_result = hotspot_result
        self.object = object
        self.tolerance = tolerance

        if isinstance(object, Protein):
            self._scored_object = self.score_protein()

        elif isinstance(object, Molecule):
            self._scored_object = self.score_molecule()

        elif isinstance(object, Cavity):
            self._scored_object = self.score_cavity()

        elif not object:
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
        :return:
        """
        # TODO: enable cavities to be generated from Protein objects
        #
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
                grid_type = interaction_pairs[feature.type]

                if feature.type == "aliphatic" or feature.type == "aromatic":
                    coordinates = Helper.cavity_centroid(feature)

                else:
                    v = feature.protein_vector
                    translate = tuple(map((h_bond_distance).__mul__, (v.x, v.y, v.z)))
                    c = feature.coordinates
                    coordinates = tuple(map(operator.add, (c.x, c.y, c.z), translate))

                if feature.atom:
                    score = self._score_atom_type(grid_type, coordinates)
                    prot.atoms[feature.atom.index].partial_charge = score
                    a = [a.index for a in prot.atoms[feature.atom.index].neighbours
                         if int(a.atomic_number) == 1]

                    if len(a) > 0:
                        for atm in a:
                            prot.atoms[atm].partial_charge = score

                else:
                    continue

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

    def _score_atom_type(self, grid_type, coordinates):
        """
        atom
        :param grid_type:
        :param coordinate:
        :param tolerance:
        :return:
        """
        if grid_type == "doneptor":
            grid_type = self._doneptor_grid(coordinates)

        return self.hotspot_result.super_grids[grid_type].value_at_coordinate(coordinates,
                                                                              tolerance=self.tolerance,
                                                                              position=False)

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

    def _doneptor_grid(self, coordinates):
        """
        An atom is scored from the grid which yields the highest value
        :param coordinates:
        :param grid_type:
        :return:
        """
        scores = [self.hotspot_result.super_grids["donor"].value_at_coordinate(coordinates,
                                                                               tolerance=self.tolerance,
                                                                               position=False),
                  self.hotspot_result.super_grids["acceptor"].value_at_coordinate(coordinates,
                                                                                  tolerance=self.tolerance,
                                                                                  position=False)
                  ]
        d = dict(zip(scores, ["donor", "acceptor"]))
        return d[max(d.keys())]

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