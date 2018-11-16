#!/usr/bin/env python

"""
More information about the fragment hotspot maps method is available from:
    Radoux, C.J. et. al., Identifying the Interactions that Determine Fragment Binding at Protein Hotspots J. Med. Chem.
    2016, 59 (9), 4314-4325
    dx.doi.org/10.1021/acs.jmedchem.5b01980
"""
from __future__ import print_function, division

import collections
import copy
import multiprocessing
import operator
import random
import sys
import tempfile
import time
from os import system, environ
from os.path import join

import nglview as nv
import numpy as np
import pkg_resources
from atomic_hotspot_calculation import AtomicHotspot, AtomicHotspotResult
from ccdc.cavity import Cavity
from ccdc.io import MoleculeWriter, MoleculeReader
from ccdc.molecule import Molecule
from ccdc.protein import Protein
from ccdc.utilities import PushDir
from grid_extension import Grid
from ipywidgets import IntSlider, interact
from hotspot_pharmacophore import PharmacophoreModel
from scipy.stats import percentileofscore
from tqdm import tqdm
from hotspot_utilities import Figures, Helper


Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


class Buriedness(object):
    """
    class to handle ghecom run
    """

    class Settings(object):
        """
        settings for ghecom run
        """

        def __init__(self):
            self.ghecom_executable = None
            self.grid_spacing = 0.5
            self.radius_min_large_sphere = 2.5
            self.radius_max_large_sphere = 9.5
            self.protein = None
            self.out_grid = None
            self.mode = "M"
            self.working_directory = tempfile.mkdtemp()
            self.in_name = join(self.working_directory, "protein.pdb")
            self.out_name = join(self.working_directory, "ghecom_out.pdb")

    def __init__(self, protein, out_grid, settings=None):
        """"""

        if settings is None:
            settings = self.Settings()

        self.settings = settings
        self.settings.protein = protein
        self.settings.out_grid = out_grid

    def calculate_buriedness(self):
        """
        runs ghecom in temporary directory

        :return: a :class: `GhecomResult` instance
        """

        with PushDir(self.settings.working_directory):

            if self.settings.protein is not None:
                with MoleculeWriter('protein.pdb') as writer:
                    writer.write(self.settings.protein)

            cmd = "{} {} -M {} -gw {} -rli {} -rlx {} -opoc {}".format(environ['GHECOM_EXE'],
                                                                       self.settings.in_name,
                                                                       self.settings.mode,
                                                                       self.settings.grid_spacing,
                                                                       self.settings.radius_min_large_sphere,
                                                                       self.settings.radius_max_large_sphere,
                                                                       self.settings.out_name)

            system(cmd)

        return BuriednessResult(self.settings)


class BuriednessResult(object):
    """
    class to store a Ghecom result
    """

    def __init__(self, settings):
        self.settings = settings
        if self.settings.out_grid:
            self.grid = self.settings.out_grid
        else:
            self.grid = self._initalise_grid(padding=1)
        self.update_grid()

    def _initalise_grid(self, padding=1):
        """
        install grid over protein to hold scores
        :param padding: value of buffer added to coordinate extremities
        :return: a `ccdc.utilities.Grid` instance
        """

        x = []
        y = []
        z = []

        for atm in self.protein.atoms:
            x.append(atm.coordinates.x)
            y.append(atm.coordinates.y)
            z.append(atm.coordinates.z)

        bl = (min(x) - padding, min(y) - padding, min(z) - padding)
        tr = (max(x) + padding, max(y) + padding, max(z) + padding)

        return Grid(origin=bl, far_corner=tr, spacing=0.5)

    def update_grid(self):
        """
        update initialised grid with ghecom values
        :return: None
        """

        lines = Helper.get_lines_from_file(self.settings.out_name)
        for line in lines:
            if line.startswith("HETATM"):
                coordinates = (float(line[31:38]), float(line[39:46]), float(line[47:54]))
                rinacc = float(line[61:66])
                i, j, k = self.grid.point_to_indices(coordinates)
                x, y, z = self.grid.nsteps
                if 0 < i < x and 0 < j < y and 0 < k < z:
                    self.grid.set_value(i, j, k, 9.5 - rinacc)


class _WeightedResult(object):
    """
    class to hold weighted grids
    """

    def __init__(self, identifier, grid):
        self.identifier = identifier
        self.grid = grid


class _SampleGrid(object):
    """
    class to handle sampled grids
    """

    def __init__(self, name, grid, atom_predicate):
        """
        attributes of SampleGrid
        :param name: str, name of probe (donor, acceptor, apolar, positive, negative)
        :param grid: a :class: `ccdc.utilities.Grid` instance
        :param atom_predicate: atom_predicate will be used to select atoms of a molecule for sampling
        """
        self.name = name
        self.grid = grid
        self.atom_predicate = atom_predicate
        self.mol = None

    @staticmethod
    def add_coordinates(coord, trans):
        """
        provides a coordinate list of atoms to be scored be scored in the SampleGrid
        :param coord: tup, (float(x), float(y), float(z)), set of atomic coordinates for "active" coordinates
        :param trans: tup, (float(x), float(y), float(z)), set of translations to translate probes to points
        above threshold
        :return: list of tup
        """

        return [coord[i] + trans[i] for i in xrange(0, len(trans))]

    def sample(self, coordinate_list, trans):
        """
        score Molecule in grids for which it has active atoms
        :param coordinate_list: list, set of coordinates of translated and rotated probe atoms to be scored
        :param trans:
        :return:
        """

        try:
            return [self.grid._grid.value(self.add_coordinates(c, trans)) for c in coordinate_list]
        except RuntimeError:
            return [0]

    def set_molecule(self, mol, polar_contribution):
        """
        set which atoms match the grid
        probes with polar atoms contribute do not contribute to apolar maps as this leads to artefacts

        :param mol:
        :param polar_contribution:
        :return:
        """
        self.mol = mol
        if self.name == 'apolar' and not polar_contribution and len([a for a in mol.atoms
                                                                     if a.is_donor or a.is_acceptor]) > 0:
            self._active_atoms = []
        elif not hasattr(self, '_active_atoms'):
            self._active_atoms = [a for a in self.mol.atoms if self.atom_predicate(a)]

    @staticmethod
    def is_donor(a):
        """
        returns true if a given atom is a donor

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom classification is "donor"
        """

        if a.is_donor and a.formal_charge == 0:
            return True
        else:
            return False

    @staticmethod
    def is_acceptor(a):
        """
        returns true if a given atom is a acceptor

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom classification is "acceptor"
        """
        if a.is_acceptor and a.formal_charge == 0:
            return True
        else:
            return False

    @staticmethod
    def is_apolar(a):
        """
        returns true if a given atom is a apolar

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom classification is "apolar"
        """
        if a.is_donor or a.is_acceptor or a.formal_charge != 0 or a.atomic_symbol == "Xe":
            return False
        else:
            return True

    @staticmethod
    def is_positive(a):
        """
        returns true if a given atom is a positively charged

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom is positively charged
        """
        return a.formal_charge > 0

    @staticmethod
    def is_negative(a):
        """
        returns true if a given atom is a negatively charged

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom is negatively charged
        """
        return a.formal_charge < 0

    @staticmethod
    def is_aromatic(a):
        """
        returns true if a given atom is aromatic

        :param a: a `ccdc.molecule.Atom` instance
        :return: bool, true if the atom is aromatic
        """
        return a.atomic_symbol == 'C' and a.is_cyclic and any(b.atom_type == 'aromatic' for b in a.bonds)


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

            values = [a.partial_charge for a in self.scored_object.heavy_atoms]
            self._score = self._geometric_mean(values=values)

        elif isinstance(object, Cavity):
            self._scored_object = self.score_cavity()

        elif not object:
            self._scored_object = self.score_hotspot()

        else:
            raise IOError("supplied object not currently supported, soz!")

    @property
    def scored_object(self):
        return self._scored_object

    @property
    def score(self):
        return self._score

    def score_protein(self):
        """
        score a protein's atoms, values stored as partial charge
        h_bond_distance between 1.5 - 2.5 A (2.0 A taken for simplicity)
        :return:
        """
        # TODO: enable cavities to be generated from Protein objects
        #
        prot = copy.copy(self.object)
        h_bond_distance = 2.0
        interaction_pairs = {"acceptor": "donor",
                             "donor": "acceptor",
                             "pi": "apolar",
                             "aliphatic":"apolar",
                             "aromatic": "apolar",
                             "apolar": "apolar",
                             "donor_acceptor": "doneptor",
                             "dummy": "dummy"}

        cavities = Helper.cavity_from_protein(self.hotspot_result.protein)
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
                else:
                    print("WARNING: no atom")

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
        # TODO: return scored cavity features, the score protein function should be enough tbh
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
            coordinates = atom.coordinates
            if atom_type == "doneptor":
                atom_type = self._doneptor_grid(atom.coordinates, grid_type=True)
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
    def _geometric_mean(values):
        '''Calculate geometric mean of scores'''
        return reduce(operator.__mul__, values, 1.0) ** (1. / len(values))

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


class HotspotResults(object):
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

    class HotspotFeature(object):
        """
        class to hold polar islands above threshold "features"
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
            self._score = self.score_feature()

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
        def score(self):
            return self._score

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

    def score(self, obj=None, return_value=False, tolerance=2):
        """
        Given a supported CCDC object, will return the object annotated with Fragment Hotspot scores
        :param obj:
        :param return_value:
        :return:
        """
        scorer = _Scorer(self, obj, tolerance)

        if return_value:
            return scorer.score, scorer.scored_object

        else:
            return scorer.scored_object

    def get_selectivity_map(self, other):
        '''
        Generate maps to highlight selectivity for a target over an off target cavity. Proteins should be aligned
        by the binding site of interest prior to calculation of Fragment Hotspot Maps. High scoring regions of a map
        represent areas of favourable interaction in the target binding site, not present in off target binding site

        :param other: a :class:`fragment_hotspots.Hotspots.HotspotResults` instance
        :return: a :class:`fragment_hotspots.Hotspots.HotspotResults` instance
        '''

        selectivity_grids = {}
        for probe in self.super_grids.keys():
            g1 = self.super_grids[probe]
            g2 = other.super_grids[probe]
            og1, og2 = self._common_grid(g1, g2)
            sele = og1 - og2
            selectivity_grids[probe] = sele
        hr = Hotspots.HotspotResults(selectivity_grids, self.protein, self.fname, None, None, self.out_dir)
        return hr

    def get_pharmacophore_model(self, identifier="id_01", cutoff=5):
        """
        method of hotspot results object(intended to be run after extracted HS) returns PharmacophoreModel
        :return:
        """
        return PharmacophoreModel.from_hotspot(self.protein, self.super_grids, identifier=identifier, cutoff=cutoff)

    def get_histogram(self, fpath="histogram.png", plot=True):
        """
        get histogram data
        :param fpath: path to output file
        :param plot:
        :return:
        """
        if plot:
            data, plt = Figures.histogram(self, plot)
            plt.savefig(fpath)
            return data, plt
        else:
            data = Figures.histogram(self, plot)
            return data

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
        :param large_cavities:
        :return:
        """
        centroid = self.location.centroid()
        cavity = [c for c in large_cavities if c.contains_point(centroid)]

        if len(cavity) == 0:
            return 0
        else:
            return 1

    def get_superstar_profile(self, feature_radius=1.5, nthreads=6, features=None, best_volume=None):
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

                ss_dict = {Helper.get_distance(feat.feature_coordinates, island.centroid()) :island
                           for island in r.grid.islands(threshold=1)
                           if Helper.get_distance(feat.feature_coordinates, island.centroid()) < 1}

                if len(ss_dict) == 0:
                    g = r.grid.copy_and_clear()

                else:
                    shortest = sorted([f[0] for f in ss_dict.items()], reverse=False)[0]
                    g = ss_dict[shortest]


                feat.superstar_results.append(AtomicHotspotResult(identifier=r.identifier,
                                                                  grid= g,
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
                        f.append(HotspotResults.HotspotFeature(probe, island))
        return f

    def _rank_features(self):
        """
        rank features based upon feature score (TO DO: modify score if required)
        :return:
        """
        feature_by_score = {feat.score: feat for feat in self.features}
        score = sorted([f[0] for f in feature_by_score.items()], reverse=True)
        for i, key in enumerate(score):
            feature_by_score[key]._rank = int(i + 1)


    def extract_pocket(self, whole_residues=False):
        '''
        Create a :class:`ccdc.Protein` containing atoms or residues that have a score

        :param whole_residues: bool, whether to include all residue atoms if only a subset have a score > 0
        :return: a :class:`ccdc.Protein` instance
        '''
        prot_scores = self.score_protein()
        pocket = self.protein.copy()
        pocket.remove_hydrogens()
        for residue in pocket.residues:
            keep_residue = False
            for atom in residue.atoms:
                # if atom.atomic_number == 1:
                #     continue
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

    def ngl_widget(self, out_dir=None):
        """
        jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000000.0
        creates ngl widget from hotspot. For use in ipython notebooks
        :param out_dir:
        :return:
        """
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


def _sample_job(q, ):

    return 0


class Hotspots(object):
    """
    A class for running Fragment Hotspot Map calculations
    """

    class _Sampler(object):
        """
        Samples one or more grids with a probe molecule
        """

        class Settings(object):
            """
            Settings for the sampler
                nrotations:                   number of rotations (keep it below 10**6)
                apolar_translation_threshold: translate probe to grid points above this threshold. Give lower values for
                                              greater sampling. Default 15
                polar_translation_threshold:  translate probe to grid points above this threshold. Give lower values for
                                              greater sampling. Default 15
                polar_contributions:          allow carbon atoms of probes with polar atoms to contribute to the apolar
                                              output map.
            """

            nrotations = 3000
            apolar_translation_threshold = 15
            polar_translation_threshold = 15
            polar_contributions = False
            return_probes = False

        def __init__(self, *grids, **kw):
            """
            Settings used to run fragment-hotspot-maps script

            :param grids: list, list of :class: `ccdc.utilities.Grid` instances
            :param kw: 'settings'
            """

            self.grids = grids
            self.grid_dic = {}
            for g in grids:
                self.grid_dic[g.name] = g
            settings = kw.get('settings')
            if settings is None:
                settings = self.Settings()

            self.settings = settings
            self.probe_grids = [_SampleGrid(g.name, g.grid.copy_and_clear(), g.atom_predicate) for g in self.grids]

        def get_priority_atom(self, molecule):
            """
            Select priority atom. Select polar atom. If multiple polar atoms, select the one furthest from the centre of
            geometry. If no polar atoms, select atom furthest from centre of geometry

            :param molecule: a :class: `ccdc.molecule.Molecule` instance
            :return: a :class: `ccdc.molecule.Molecule` instance, str atom type
            """
            c = molecule.centre_of_geometry()
            polar_atoms = [a for a in molecule.atoms if a.is_donor or a.is_acceptor]
            atom_by_distance = {}
            if len(polar_atoms) > 0:
                for a in polar_atoms:
                    d = Helper.get_distance(c, a.coordinates)
                    atom_by_distance[d] = a
            else:
                for a in molecule.atoms:
                    d = Helper.get_distance(c, a.coordinates)
                    atom_by_distance[d] = a

            greatest_distance = sorted(atom_by_distance.keys())[0]
            priority_atom = atom_by_distance[greatest_distance]

            pa_type = None
            if priority_atom.formal_charge != 0:
                if priority_atom.formal_charge < 0:
                    pa_type = "negative"
                elif priority_atom.formal_charge > 0:
                    pa_type = "positive"
            else:
                if priority_atom.is_acceptor:
                    pa_type = "acceptor"
                elif priority_atom.is_donor:
                    pa_type = "donor"
                else:
                    pa_type = "apolar"

            return priority_atom, pa_type

        def get_translation_points(self, priority_atom_type):
            """
            returns a list of coordinates that are greater than the threshold, that the probe will be translated to

            :param priority_atom_type: str, atomic interaction type
            :return: list, list of :class: `ccdc.molecule.Molecule` instances
            """
            translate_probe = []
            wg = self.grid_dic[priority_atom_type]

            if priority_atom_type == 'apolar':
                translation_threshold = self.settings.apolar_translation_threshold
            else:
                translation_threshold = self.settings.polar_translation_threshold
            hs = wg.grid.islands(translation_threshold)
            for g in hs:
                nx, ny, nz = g.nsteps

                maxima = [g.indices_to_point(i, j, k)
                          for i in range(nx) for j in range(ny) for k in range(nz)
                          if g.value(i, j, k) >= translation_threshold]

                translate_probe = translate_probe + maxima
            return translate_probe

        def generate_rand_quaternions(self):
            """
            Returns a list of random quaternions. Length matches settings.nrotations

            :return: tup, (a,b,c,d)
            """
            quaternions = []
            i = 1
            if self.settings.nrotations > 1:
                while i <= self.settings.nrotations:
                    r1 = random.uniform(-1, 1)
                    r2 = random.uniform(-1, 1)
                    s1 = r1 * r1 + r2 * r2
                    if s1 < 1:
                        r3 = random.uniform(-1, 1)
                        r4 = random.uniform(-1, 1)

                        s2 = r3 * r3 + r4 * r4
                        if s2 < 1:
                            q = (r1, r2, r3 * (np.sqrt((1 - s1) / s2)), r4 * (np.sqrt((1 - s1) / s2)))
                            quaternions.append(q)

                            i += 1

            return quaternions

        @staticmethod
        def score(values):
            """
            Calculate geometric mean of scores

            :param values: float, scores of atoms in probe
            :return: float, geometric mean of probe atom scores
            """

            return reduce(operator.__mul__, values, 1.0) ** (1. / len(values))

        def sample_pose(self, trans, active_atoms_dic, probe):
            """
            Return a pose score (as defined by the score(self,dic) function) and a dictionary of atom:scores

            :param trans: list of translations
            :param active_atoms_dic: dict {"interaction type": "atoms to be scored"}
            :param probe: str, interaction_type
            :return:
            """
            if probe == "negative" or probe == "positive":
                atoms = [g.mol.atoms for g in self.grids if g.name == "positive"][0]
                weight = int(6 / len([a for a in atoms if str(a.atomic_symbol) == "C"]))

                apolar_values = [g.sample(active_atoms_dic[g.name], trans) for g in self.grids if
                                 len(active_atoms_dic[g.name]) > 0 and g.name == "apolar"] * weight
                charged_values = [g.sample(active_atoms_dic[g.name], trans) for g in self.grids if
                                  len(active_atoms_dic[g.name]) > 0 and g.name != "apolar"]
                values = apolar_values + charged_values

            else:
                values = [g.sample(active_atoms_dic[g.name], trans) for g in self.grids if
                          len(active_atoms_dic[g.name]) > 0]
            scores = [item for sublist in values for item in sublist]
            return self.score(scores)

        def update_out_grids(self, score, active_coordinates_dic, trans):
            """
            For active atoms for a given grid, set closest grid point value to score, unless already set to a higher
            value

            :param score: float, score of a given probe
            :param active_coordinates_dic:
            :param trans: list of translations
            :return:
            """

            for pg in self.probe_grids:
                actives = active_coordinates_dic[pg.name]
                if len(actives) == 0:
                    continue
                for a in actives:
                    i, j, k = pg.grid.point_to_indices(pg.add_coordinates(a, trans))
                    pg.grid.set_value(i, j, k, max(score, pg.grid.value(i, j, k)))

        def get_active_coordinates(self):
            """
            Returns a dictionary of {grid_name:[Coordinates]}

            :return:
            """
            active_coords_dic = {
                g.name: [a.coordinates for a in g._active_atoms]
                for g in self.grids
                }

            return active_coords_dic

        def sample(self, molecule, probe):
            """
            Sample the grids according to the settings

            :param molecule:
            :param probe: str, interaction type, (donor, acceptor, negative, positive, apolar)
            :return:
            """
            high_scoring_probes = {}
            priority_atom, priority_atom_type = self.get_priority_atom(molecule)
            translate_points = self.get_translation_points(priority_atom_type)
            molecule.remove_hydrogens()
            quaternions = self.generate_rand_quaternions()

            print("\n    nRotations:", len(quaternions), "nTranslations:" , len(translate_points), "probename:", probe)

            for g in self.grids:
                g.set_molecule(molecule, True)

            for g in self.probe_grids:
                g.set_molecule(molecule, self.settings.polar_contributions)

            for q in tqdm(quaternions):
                molecule.apply_quaternion(q)
                priority_atom_coordinates = priority_atom.coordinates
                active_coordinates_dic = self.get_active_coordinates()

                for priority_atom_point in translate_points:

                    translation = [priority_atom_point[i] - priority_atom_coordinates[i]
                                   for i in xrange(0, len(priority_atom_coordinates))]

                    score = self.sample_pose(translation, active_coordinates_dic, probe)
                    self.update_out_grids(score, active_coordinates_dic, translation)

                    if self.settings.return_probes:
                        if score < 5:
                            continue
                        if score > 14:
                            m = molecule.copy()
                            m.translate(translation)
                            m.identifier = "{}".format(score)

                            try:
                                high_scoring_probes[score].append(m)
                            except KeyError:
                                high_scoring_probes[score] = [m]

            if self.settings.return_probes:
                sampled_probes = []
                for key in sorted(high_scoring_probes.iterkeys(), reverse=True):
                    sampled_probes.extend(high_scoring_probes[key])

                if len(sampled_probes) > 10000:
                    return sampled_probes[:10000]
                else:
                    return sampled_probes

    def __init__(self, settings=None):
        self.out_grids = {}
        self.super_grids = {}

        # self.superstar_grids = None
        # self.weighted_grids = None
        # self.sampled_probes = {}
        #
        # self.protein = None
        # self.fname = None
        # self.probe_size = None
        # self.charged_probes = None
        #
        # #self.out_dir = None
        # self.wrk_dir = None

        if settings is None:
            self.sampler_settings = self._Sampler.Settings()
        else:
            self.sampler_settings = settings

    @property
    def protein(self):
        return self._protein

    @protein.setter
    def protein(self, prot):
        if isinstance(prot, Protein):
            self._protein = prot
        else:
            raise TypeError("`ccdc.protein.Protein` must be supplied. Hint: Use Protein.from_file()")

    @property
    def charged_probes(self):
        return self._charged_probes

    @charged_probes.setter
    def charged_probes(self, option):
        if type(option) is bool:
            self._charged_probes = option
        else:
            raise TypeError("Expecting a bool, got {} instead".format(type(option)))

    @property
    def probe_size(self):
        return self._probe_size

    @probe_size.setter
    def probe_size(self, size):
        if size in range(3,8):
            self._probe_size = size
        else:
            raise ValueError("Probe size must be an integer between 3-7")

    @property
    def buriedness_method(self):
        return self._buriedness_method

    @buriedness_method.setter
    def buriedness_method(self, method):
        method = method.lower()
        if method == 'ghecom':
            if sys.platform == 'linux' or sys.platform == 'linux2':
                if 'GHECOM_EXE' in environ:
                    self._buriedness_method = method
                else:
                    raise EnvironmentError("Must set Ghecom environment variable")
            else:
                raise OSError('Ghecom is only supported on linux')

        elif method == 'ligsite':
            if sys.platform == 'linux' or sys.platfrom == 'linux2':
                print("RECOMMENDATION: you have chosen LIGSITE as buriedness method, ghecom is recommended")
            self._buriedness_method = method

        else:
            raise ValueError("Buriedness method must be 'ghecom' (default) or 'ligsite")

    @property
    def cavities(self):
        return self._cavities

    @cavities.setter
    def cavities(self, obj):
        if obj is not None:
            if isinstance(obj, list) or isinstance(obj, tuple):
                if isinstance(obj, Coordinates):
                    self._cavities = obj
                elif isinstance(obj, Molecule):
                    self._cavities = [m.centre_of_geometry() for m in obj]
                elif isinstance(obj, Cavity):
                    self._cavities = [Helper.cavity_centroid(c) for c in obj]
                else:
                    self._cavities = None

            elif isinstance(obj, Coordinates):
                self._cavities = [obj]
            elif isinstance(obj, Molecule):
                self._cavities = [obj.centre_of_geometry()]
            elif isinstance(obj, Cavity):
                self._cavities = [Helper.cavity_centroid(obj)]

            else:
                raise TypeError("Type unsupported. Hint: Cavity can be list, Coordinate, Molecule or Cavity")

    @property
    def nprocesses(self):
        return self._nprocesses

    @nprocesses.setter
    def nprocesses(self, num):
        num = int(num)
        if num in range(0, int(multiprocessing.cpu_count())):
            self._nprocesses = num
        else:
            raise OSError("CPU count = {}".format(multiprocessing.cpu_count()))

    @property
    def sampler_settings(self):
        return self._sampler_settings

    @sampler_settings.setter
    def sampler_settings(self, settings):
        if isinstance(settings, self._Sampler.Settings):
            self._sampler_settings = settings
        else:
            self._sampler_settings = None

    def _get_weighted_maps(self):
        """
        weight superstar output by burriedness
        :return: a list of :class: `WeightedResult` instances
        """
        results = []
        for s in self.superstar_grids:
            g, b = Grid.common_grid([s.grid, s.buriedness], padding=1)
            weighted_grid = g * b
            results.append(_WeightedResult(s.identifier, weighted_grid))

        return results

    def _get_out_maps(self, probe, grid_dict, return_probes=False):
        """
        Function to organise sampling of weighted superstar maps by molecular probes

        :param probe:
        :param grid_dict:
        :return:
        """

        donor_grid = _SampleGrid('donor', grid_dict['donor'], _SampleGrid.is_donor)
        acceptor_grid = _SampleGrid('acceptor', grid_dict['acceptor'], _SampleGrid.is_acceptor)
        apolar_grid = _SampleGrid('apolar', grid_dict['apolar'], _SampleGrid.is_apolar)

        if self.charged_probes:
            negative_grid = _SampleGrid('negative', grid_dict['negative'], _SampleGrid.is_negative)
            positive_grid = _SampleGrid('positive', grid_dict['positive'], _SampleGrid.is_positive)

        kw = {'settings': self.sampler_settings}
        if self.charged_probes:
            self.sampler = self._Sampler(apolar_grid, donor_grid, acceptor_grid, negative_grid, positive_grid, **kw)
        else:
            self.sampler = self._Sampler(apolar_grid, donor_grid, acceptor_grid, **kw)

        probe_path = pkg_resources.resource_filename('hotspots', 'probes/')

        if self.charged_probes:
            if probe == "negative" or probe == "positive":
                mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, "test")))[0]
            else:
                mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, self.probe_size)))[0]
        else:
            mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, self.probe_size)))[0]

        probes = self.sampler.sample(mol, probe=probe, return_probes=return_probes)

        for pg in self.sampler.probe_grids:
            if pg.name.lower() == probe:
                try:
                    self.out_grids[pg.name].append(pg.grid)
                except KeyError:
                    self.out_grids[pg.name] = [pg.grid]

        if return_probes:
            return probes

    def _calc_hotspots(self, return_probes=False):
        """
        Function for overall organisation of hotspot calculation
        :return:
        """
        print("Start atomic hotspot detection")
        a = AtomicHotspot()
        a.settings.atomic_probes = {"apolar" : "AROMATIC CH CARBON",
                                    "donor" : "UNCHARGED NH NITROGEN",
                                    "acceptor" : "CARBONYL OXYGEN"}
        if self.charged_probes:
            a.settings.atomic_probes = {"negative": "CARBOXYLATE OXYGEN", "positive": "CHARGED NH NITROGEN"}

        probe_types = a.settings.atomic_probes.keys()
        self.superstar_grids = a.calculate(protein=self.protein,
                                           nthreads=self.nprocesses,
                                           cavity_origins=self.cavities)

        print("Atomic hotspot detection complete\n")

        print("Start buriedness calcualtion")
        if self.buriedness_method == 'ghecom':
            print("    method: Ghecom")
            out_grid = self.superstar_grids[0].buriedness.copy_and_clear()
            b = Buriedness(protein=self.protein,
                           out_grid=out_grid)
            self.buriedness = b.calculate_buriedness().grid
        else:
            print("    method: LIGSITE")
            self.buriedness = Grid.get_single_grid(grd_dict={s.identifier: s.buriedness for s in self.superstar_grids},
                                                   mask=False)

        self.weighted_grids = self._get_weighted_maps()

        print("Buriedness calcualtion complete\n")

        print("Start sampling")
        grid_dict = {w.identifier: w.grid for w in self.weighted_grids}

        for probe in probe_types:
            if return_probes:
                self.sampled_probes.update(probe, self._get_out_maps(probe, grid_dict))

            else:
                self._get_out_maps(probe, grid_dict)

        print("Sampling complete\n")

    def from_protein(self, protein, charged_probes=False, probe_size=7, buriedness_method= 'ghecom',
                     cavities=None, nprocesses=None, sampler_settings=None):
        """

        :param protein: a :class:`ccdc.protein.Protein` instance
        :param charged_probes: bool, if True include positive and negative probes
        :param probe_size: int, size of probe in number of heavy atoms (3-8 atoms)
        :param buriedness_method: str, either 'ghecom' or 'ligsite'
        :param cavities: Coordinate or `ccdc.cavity.Cavity` or `ccdc.molecule.Molecule` or list,
        algorithm run on parsed cavity
        :param nprocesses: int, number of CPU's used
        :param sampler_settings: a `hotspots.Hotspot._Sampler.Settings`, holds the sampler settings
        :return:
        """

        start = time.time()
        self.protein = protein
        self.charged_probes = charged_probes
        self.probe_size = probe_size
        self.buriedness_method = buriedness_method
        self.cavities = cavities
        self.nprocesses = nprocesses
        self.sampler_settings = sampler_settings

        self._calc_hotspots()     # return probes = False by default
        self.super_grids = {p: g[0] for p, g in self.out_grids.items()}

        print("Runtime = {}seconds".format(time.time() - start))

        return HotspotResults(super_grids=self.super_grids,
                              protein=self.protein,
                              buriedness=self.buriedness)
