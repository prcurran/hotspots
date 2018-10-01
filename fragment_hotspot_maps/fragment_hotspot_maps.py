#!/usr/bin/env python

"""
More information about the fragment hotspot maps method is available from:
    Radoux, C.J. et. al., Identifying the Interactions that Determine Fragment Binding at Protein Hotspots J. Med. Chem.
    2016, 59 (9), 4314-4325
    dx.doi.org/10.1021/acs.jmedchem.5b01980
"""
from __future__ import print_function, division

import argparse
import math
import operator

from os.path import join, dirname, exists, splitext
from os import environ, name, getcwd, mkdir, chdir, system, path, listdir
import sys
import glob
import random
import subprocess
import zipfile
import shutil
import tempfile
import json

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from ccdc.io import csd_directory, MoleculeWriter, MoleculeReader
from ccdc.protein import Protein
from ccdc.utilities import _test_output_dir, PushDir
from enhanced_grid import Grid
from ccdc import molecule

from scipy import optimize, ndimage
from skimage import feature
import pkg_resources
from template_strings import pymol_template, colourmap, superstar_ins, extracted_hotspot_template, pymol_zip_template, \
    crossminer_header
from extraction import BuildLocation, BuildFeature, HotspotBuilder
from pharmacophore import PharmacophoreModel
import pandas as pd

from concurrent import futures
from pprint import pprint

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import AllChem
    from matplotlib.colors import LinearSegmentedColormap
except ImportError:
    print("""ImportError: rdkit(optional)\nInfo: This module is required for producing 2D schematic maps""")
if name == 'nt':
    pass


class _HotspotsHelper(object):
    """
    class providing utility functions
    """

    @staticmethod
    def _indices_to_point(i, j, k, g):
        """
        Return x,y,z coordinate for a given grid index

        :param i: int, indice (x-axis)
        :param j: int, indice (y-axis)
        :param k: int, indice (z-axis)
        :param g: a :class: `ccdc.utilities.Grid` instance
        :return: float(x), float(y), float(z)
        """

        ox, oy, oz, = g.bounding_box[0]
        gs = 0.5
        return ox + float(i) * gs, oy + float(j) * gs, oz + gs * float(k)

    @staticmethod
    def _point_to_indices(p, g):
        """
        Return the nearest grid index for a given point

        :param p: tup, (float(x), float(y), float(z)), a coordinate on a 3D grid
        :param g: a :class: `ccdc.utilities.Grid` instance
        :return: int(x), int(y), int(z)
        """

        gs = 0.5
        rx, ry, rz = [round(i / gs) for i in p]
        ox, oy, oz = [round(i / gs) for i in g.bounding_box[0]]
        return int(rx - ox), int(ry - oy), int(rz - oz)

    @staticmethod
    def _get_distance(coords1, coords2):
        """
        given two coordinates, calculates the distance

        :param coords1: tup, (float(x), float(y), float(z), coordinates of point 1
        :param coords2: tup, (float(x), float(y), float(z), coordinates of point 2
        :return: float, distance
        """
        xd = coords1[0] - coords2[0]
        yd = coords1[1] - coords2[1]
        zd = coords1[2] - coords2[2]
        d = math.sqrt(xd ** 2 + yd ** 2 + zd ** 2)
        return d

    @staticmethod
    def _copy_and_clear(grid):
        """
        make a new empty grid

        :param grid: a :class: `ccdc.utilities.Grid` instance
        :return: an empty :class: `ccdc.utilities.Grid` instance
        """

        g = grid.copy()
        g *= 0
        return g

    def _common_grid(self, g1, g2, padding=1):
        """
        finds a common frame of reference for two grids

        :param g1: a :class: `ccdc.utilities.Grid` instance
        :param g2: a :class: `ccdc.utilities.Grid` instance
        :param padding: int, additional grid point in x, y, z directions
        :return:
        """

        sg = Grid.super_grid(padding, g1, g2)
        out_g = self._copy_and_clear(sg)
        out1 = Grid.super_grid(padding, g1, out_g)
        out2 = Grid.super_grid(padding, g2, out_g)
        return out1, out2

    @staticmethod
    def contains_point(point, island, tolerance=0, threshold=0):

        mini = island.bounding_box[0]
        maxi = island.bounding_box[1]
        if island.value_at_point(point) > threshold:
            return all([mini.x - tolerance < point[0] < maxi.x + tolerance,
                        mini.y - tolerance < point[1] < maxi.y + tolerance,
                        mini.z - tolerance < point[2] < maxi.z + tolerance])
        else:
            return False

    @staticmethod
    def grid_score(g, threshold=0, percentile=75):
        """
        take a group and return average score of points above threshold
        :param g:
        :param threshold
        :param percentile

        :return:
        """

        nx, ny, nz = g.nsteps
        values = [g.value(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)
                  if g.value(i, j, k) > threshold]

        if len(values) == 0:
            return 0
        else:
            return np.percentile(values, percentile)

    @staticmethod
    def _get_lines_from_file(fname):
        """
        fetch lines from ghecom output

        :return: str, lines from output file
        """

        f = open(fname)
        lines = f.readlines()
        f.close()
        for i in range(0, len(lines)):
            lines[i] = lines[i].strip()
        return lines


class _RunSuperstar(object):
    """
    class to handle SuperStar run
    """

    class Settings(object):
        """
        setting for Superstar run
        """

        def __init__(self):
            self.jobname = None
            self.probename = None
            self.moleculefile = None
            self.cavity_origin = None

            # self.occulsionthreshold = 5
            self.mapbackgroundvalue = 1
            self.boxborder = 10
            self.minpropensity = 1
            self.superstar_executable = None
            self.superstar_env = None
            self.working_directory = None

    def __init__(self, **kw):
        settings = kw.get('settings')
        if settings is None:
            settings = self.Settings()
        self.settings = settings

        main_dir = environ.get('MAINDIR')
        if main_dir:
            if sys.platform == 'win32':
                self.settings.superstar_executable = 'superstar_app.exe'
            else:
                self.settings.superstar_executable = ' '.join([join(environ['MAINDIR'], 'run.sh'), 'superstar_app.x'])
            self.settings.superstar_env = dict()
        else:
            if sys.platform == 'win32':
                base = dirname(csd_directory())
                merc = glob.glob(join(base, 'mercury*'))
                if len(merc):
                    merc = merc[0]
                #self.settings.superstar_executable = join(merc, 'superstar_app.exe')
                self.settings.superstar_executable = join(merc, 'superstar*')
                self.settings.superstar_env = dict(
                    SUPERSTAR_ISODIR=str(join(dirname(csd_directory()), 'isostar_files', 'istr')),
                    SUPERSTAR_ROOT=str(join(dirname(csd_directory()), "Mercury"))
                )

            elif sys.platform == 'darwin':
                print("OS X not supported")

            else:
                self.settings.superstar_executable = join(dirname(csd_directory()), 'bin', 'superstar')

                self.settings.superstar_env = dict(
                    SUPERSTAR_ISODIR=str(join(dirname(csd_directory()), 'isostar_files', 'istr')),
                    SUPERSTAR_ROOT=str(dirname(csd_directory()))
                )
                # self.settings.working_directory = _test_output_dir()
                # print(self.settings.working_directory)

    def _append_cavity_info(self):
        """
        updates ins file with any cavity information

        :return: None
        """

        if self.settings.cavity_origin is not None:
            pnt = self.settings.cavity_origin
            extension = '\nCAVITY_ORIGIN {} {} {}'.format(pnt[0], pnt[1], pnt[2])
        else:
            extension = '\nSUBSTRUCTURE ALL'
        self.ins += extension

    def _get_inputs(self, out_dir):
        """
        assembles the ins files, uses a template string from template_strings.py

        :param out_dir: str, output directory
        :return: None
        """

        self.ins = superstar_ins(self.settings)
        self._append_cavity_info()
        out = join(out_dir, "ins")
        try:
            mkdir(out)
        except OSError:
            pass

        self.fname = join(out, "superstar_{}.ins".format(self.settings.jobname.split(".")[0]))
        with open(self.fname, "w") as w:
            w.write(self.ins)

    def run_superstar(self, prot, out_dir):
        """
        calls SuperStar as command-line subprocess

        :param prot: a :class:`ccdc.protein.Protein` instance
        :param out_dir: str, output directory
        :return:
        """

        with PushDir(self.settings.working_directory):
            if prot is not None:
                with MoleculeWriter(join('protein.pdb')) as writer:
                    writer.write(prot)
            self._get_inputs(out_dir)
            env = environ.copy()
            env.update(self.settings.superstar_env)
            cmd = self.settings.superstar_executable + ' ' + self.fname
            subprocess.call(cmd, shell=sys.platform != 'win32', env=env)
        return _SuperstarResult(self.settings)


class _SuperstarResult(object):
    """
    class to store a SuperStar result
    """

    def __init__(self, settings):
        self.settings = settings
        self.identifier = settings.jobname.split(".")[0]

        grid_path = join(self.settings.working_directory, self.identifier + ".ins.acnt")
        if exists(grid_path):
            self.grid = Grid.from_file(grid_path)
        else:
            raise AttributeError('{} superstar grid could not be found'.format(self.identifier))

        ligsite_path = join(self.settings.working_directory, self.identifier + ".ins.ligsite.acnt")
        if exists(grid_path):
            l = Grid.from_file(ligsite_path)
            self.ligsite = self.correct_ligsite(self.grid, l)
        else:
            raise AttributeError('{} ligsite grid could not be found'.format(self.identifier))

    @staticmethod
    def correct_ligsite(g, l):
        """
        Grid points where ligsite has a score of 0 (i.e. a clash) and SuperStar has a favourable score, set the ligsite
        grid point to its maximum neighbour

        :param g:
        :param l:
        :return:
        """

        mask = ((l < 1) & (g > 2))
        lc = l.copy()
        lc = lc.max_value_of_neighbours()
        correction = mask * lc
        correction = correction.mean_value_of_neighbours()
        corrected = l + correction
        return corrected


class _RunGhecom(object):
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
            self.prot = None
            self.out_grid = None
            self.mode = "M"
            self.working_directory = _test_output_dir()
            self.in_name = join(self.working_directory, "protein.pdb")
            self.out_name = join(self.working_directory, "ghecom_out.pdb")

    def __init__(self, **kw):
        """

        :param kw: settings for ghecom run
        """

        settings = kw.get('settings')
        if settings is None:
            settings = self.Settings()

        self.settings = settings
        self.settings.prot = None
        self.settings.out_grid = None

    def run_ghecom(self):
        """
        runs ghecom in temporary directory

        :return: a :class: `GhecomResult` instance
        """

        with PushDir(self.settings.working_directory):
            if self.settings.prot is not None:
                with MoleculeWriter('protein.pdb') as writer:
                    writer.write(self.settings.prot)

            cmd = "{}/ghecom {} -M {} -gw {} -rli {} -rlx {} -opoc {}".format(self.settings.ghecom_executable,
                                                                              self.settings.in_name,
                                                                              self.settings.mode,
                                                                              self.settings.grid_spacing,
                                                                              self.settings.radius_min_large_sphere,
                                                                              self.settings.radius_max_large_sphere,
                                                                              self.settings.out_name)
            #
            system(cmd)

        return _GhecomResult(self.settings)


class _GhecomResult(_HotspotsHelper):
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

        for atm in self.prot.atoms:
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

        lines = self._get_lines_from_file(self.settings.out_name)
        for line in lines:
            if line.startswith("HETATM"):
                coordinates = (float(line[31:38]), float(line[39:46]), float(line[47:54]))
                rinacc = float(line[61:66])
                i, j, k = self._point_to_indices(coordinates, self.grid)
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

class _SampleGrid(_HotspotsHelper):
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


class HotspotResults(_HotspotsHelper):
    """
    A Hotspot_results object is returned at the end of a Hotspots calculation. It contains functions for accessing
    and using the results.
    """

    def __init__(self, grid_dict, protein, fname, sampled_probes, buriedness, out_dir):
        try:
            self.super_grids = grid_dict
            for probe, g in grid_dict.items():
                b = g.bounding_box
        except:
            raise TypeError("Not a valid Grid")

        self.prot = protein
        self.fname = fname
        self.buriedness = buriedness
        self.out_dir = out_dir
        self.features_by_score = {}
        self.donor_scores = None
        self.acceptor_scores = None
        self.apolar_scores = None
        self.archive_name = None
        self.archive_loc = None
        # self.sampled_probes = self.filter_by_score(sampled_probes)
        self.sampled_probes = sampled_probes
        self.pharmacophore = None

    def get_pharmacophore_model(self, identifier="id_01", cutoff=5):
        """
        method of hotspot results object(intended to be run after extracted HS) returns PharmacophoreModel
        :return:
        """
        return PharmacophoreModel.from_hotspot(self.prot, self.super_grids, identifier=identifier, cutoff=5)

    def filter_by_score(self, sampled_probes, score=16):
        """
        filter sampled probes by score (reduce searching)
        :param sampled_probes: dict
        :param score:
        :return:
        """
        if sampled_probes is not None:
            return {probe: [m for m in mols if float(m.identifier) > score] for probe, mols in
                    sampled_probes.items()}
        else:
            return None

    def _combine(self):
        """
        combines multiple grid objects in a single grid
        :return: a :class: `ccdc.utilities.Grid` instance
        """

        sg = Grid.super_grid(0, *self.super_grids.values())
        out_g = self._copy_and_clear(sg)
        return {probe: Grid.super_grid(1, g, out_g) for probe, g in self.super_grids.items()}

    def _minimal_grid(self):
        """
        takes a results object and produces minimal common grid. (reduces memory required to store grids)

        :return: a :class: `ccdc.utilities.Grid` instance
        """

        for probe, g in self.super_grids.items():
            self.super_grids[probe] = Grid.super_grid(1, *g.islands(threshold=1))

        self.super_grids = self._combine()

    def _histogram_info(self, data, key, n):
        """
        data for the matplotlib figure

        :param data:
        :param key:
        :param n:
        :return:
        """

        colour_dict = {"acceptor": "r", "donor": "b", "apolar": "y"}
        hist, bin_edges = np.histogram(data[key], bins=range(0, 40), normed=True)
        plt.bar(bin_edges[:-1], hist, width=1, color=colour_dict[key])
        plt.xlim(min(bin_edges), max(bin_edges))
        plt.ylim(0, 0.35)
        plt.yticks([])
        if n == 0:
            plt.title("Fragment-hotspot Maps")
        if n < 2:
            plt.xticks([])
        if n == 2:
            plt.xlabel("Fragment hotspot score")
        if n == 1:
            plt.ylabel("Frequency")

    def _generate_histogram(self, data):
        """
        initialise the matplotlib figure to output histograms

        :param data:
        :return:
        """

        plt.figure(1)
        for n, key in enumerate(data.keys()):
            plt.subplot(3, 1, (n + 1))
            self._histogram_info(data, key, n)
        plt.savefig("hotspot_histogram.png")
        plt.close()

    def hotspot_histograms(self, ligand=None):
        """
        Outputs histograms of hotspot score.

        :param ligand:
        :return: None
        """

        self.data = {}
        if ligand:
            self.bs_centroid = ligand.centre_of_geometry()
            for g in self.super_grids.keys():
                grd = self.super_grids[g]
                fragment_centroid = self._point_to_indices(self.bs_centroid, grd)
                n = 8
                rx = range(fragment_centroid[0] - n, fragment_centroid[0] + n)
                ry = range(fragment_centroid[1] - n, fragment_centroid[1] + n)
                rz = range(fragment_centroid[2] - n, fragment_centroid[2] + n)

                self.data[g] = np.array(
                    [grd.value(i, j, k) for i in rx for j in ry for k in rz if grd.value(i, j, k) != 0])

        else:
            for g in self.super_grids.keys():
                grd = self.super_grids[g]
                nx, ny, nz = grd.nsteps
                self.data[g] = np.array(
                    [grd.value(i, j, k) for i in xrange(0, nx) for j in xrange(0, ny) for k in xrange(0, nz) if
                     grd.value(i, j, k) != 0])

        self._generate_histogram(self.data)

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
        hr = Hotspots.HotspotResults(selectivity_grids, self.prot, self.fname, None, None, self.out_dir)
        return hr

    def _get_near_score(self, coordinates, atom_type, tolerance):
        '''Searches nearby grid points and returns the maximum score'''

        i, j, k = self._point_to_indices(coordinates, self.super_grids[atom_type])
        nx, ny, nz = self.super_grids[atom_type].nsteps
        if nx - i < tolerance + 1 or ny - j < tolerance + 1 or nz - k < tolerance + 1:
            return 0, "outside"
        if i < tolerance or j < tolerance or k < tolerance:
            return 0, "outside"

        g = self.super_grids[atom_type]

        scores = {g.value(i + di, j + dj, k + dk): (i + di, j + dj, k + dk) for di in
                  range(-tolerance, +tolerance + 1)
                  for dj in range(-tolerance, +tolerance + 1) for dk in range(-tolerance, +tolerance + 1)}

        score = sorted(scores.keys(), reverse=True)[0]
        if score < 0.01:
            return 0, 0
        i, j, k = scores[score]
        return score, self._indices_to_point(i, j, k, g)

    def _get_atom_type(self, atom):
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

    def _update_score_dic(self, atoms_by_score, atom, residue, score, atom_type):
        try:
            atoms_by_score[score].append((atom, residue, atom_type))
        except KeyError:
            atoms_by_score[score] = [(atom, residue, atom_type)]

    def _score_protein_atoms(self):
        '''Assigns a score to each protein atom. Donor scores are assigned to polar hydrogens, rather than the heavy
        atom. Returns a dictionary of {score:[atoms]}'''
        interaction_partner_dic = {'donor': 'acceptor', 'acceptor': 'donor', 'apolar': 'apolar'}
        atoms_by_score = {}
        residue = "Blah"
        for residue in self.prot.residues:
            atom_scores = []
            for atom in residue.atoms:
                if atom.atomic_number == 1:
                    continue
                atom_type = self._get_atom_type(atom)
                if atom_type == 'donor':
                    for n in atom.neighbours:
                        if n.atomic_number == 1:
                            score = self._get_near_score(atom.coordinates, 'acceptor', tolerance=5)
                            self._update_score_dic(atoms_by_score, n, residue, score, atom_type)
                elif atom_type == 'doneptor':
                    score = self._get_near_score(atom.coordinates, 'donor', tolerance=5)
                    self._update_score_dic(atoms_by_score, atom, residue, score, 'acceptor')
                    for n in atom.neighbours:
                        if n.atomic_number == 1:
                            score = self._get_near_score(atom.coordinates, 'acceptor', tolerance=5)
                            self._update_score_dic(atoms_by_score, n, residue, score, 'donor')
                else:
                    score = self._get_near_score(atom.coordinates, interaction_partner_dic[atom_type],
                                                 tolerance=5)
                    self._update_score_dic(atoms_by_score, atom, residue, score, atom_type)

        return atoms_by_score

    def score_protein(self):
        '''
        Assigns a score to each protein atom. Donor scores are assigned to polar hydrogens, rather than the heavy
        atom.
        atom_id = "{0}/{1}/{2}".format(residue.chain_identifier, residue.identifier.split(':')[1][3:],atom.label)

        :return: dict of {atom_id: score}

        '''

        atom_dic = {}
        donor_scores = []
        acceptor_scores = []
        apolar_scores = []
        for residue in self.prot.residues:
            for atom in residue.atoms:
                if atom.atomic_number == 1:
                    continue

                donor, donor_coord = self._get_near_score(atom.coordinates, 'acceptor', tolerance=4)
                acceptor, acceptor_coord = self._get_near_score(atom.coordinates, 'donor', tolerance=4)
                apolar, apolar_coord = self._get_near_score(atom.coordinates, 'apolar', tolerance=4)
                id = "{0}/{1}/{2}".format(residue.chain_identifier, residue.identifier.split(':')[1][3:],
                                          atom.label)
                atom_dic[id] = {
                    'donor': donor,
                    'acceptor': acceptor,
                    'apolar': apolar,
                    'residue': residue,
                    'donor_coord': donor_coord,
                    'acceptor_coord': acceptor_coord,
                    'apolar_coord': apolar_coord,
                }

                if donor > 0:
                    donor_scores.append(donor)
                if acceptor > 0:
                    acceptor_scores.append(acceptor)
                if apolar > 0:
                    apolar_scores.append(apolar)

        self.donor_scores = donor_scores
        self.acceptor_scores = acceptor_scores
        self.apolar_scores = apolar_scores

        return atom_dic

    def _get_percentiles(self, percentile=90):
        if self.apolar_scores is None:
            self.score_protein()

        try:
            donor_percentile = np.percentile(self.donor_scores, percentile)
        except IndexError:
            donor_percentile = 14

        try:
            acceptor_percentile = np.percentile(self.acceptor_scores, percentile)
        except IndexError:
            acceptor_percentile = 14

        try:
            apolar_percentile = np.percentile(self.apolar_scores, percentile)
        except:
            apolar_percentile = 14

        return {'donor': donor_percentile,
                'acceptor': acceptor_percentile,
                'apolar': apolar_percentile}

    def _get_feat_atom(self, coords):
        ats = [a for a in self.prot.atoms if
               round(a.coordinates[0], 2) == round(coords[0], 2) and round(a.coordinates[1], 2) == round(coords[1],
                                                                                                         2)]
        # and round(a.coordinates[2],2) == coords[2]]
        return ats[0]

    def _get_ideal_coordinates(self, vector, feature_coords, feat_distance):
        feat_x = feature_coords.x + (vector[0] * feat_distance)
        feat_y = feature_coords.y + (vector[1] * feat_distance)
        feat_z = feature_coords.z + (vector[2] * feat_distance)

        return (feat_x, feat_y, feat_z)

    def _score_cavity_features(self, cav):
        '''
        Assign a score to each cavity feature, based on the feature's ideal interaction point

        :param cav: a :class:`ccdc.Cavity.cavity` instance
        :return: a dictionary {score:[list of features], where each item in list of features is a tuple of (atom, feature, atom_type)
        '''

        '''Assigns a score to each cavity feature. Donor scores are assigned to polar hydrogens, rather than the heavy
            atom. Returns a dictionary of {score:[features]}'''

        for feat in cav.features:
            atom_type_dict = {'pi': 'apolar', 'aromatic': 'apolar', 'aliphatic': 'apolar', 'donor': 'acceptor',
                              'acceptor': 'donor'}
            feat_coords = feat.coordinates
            atom = self._get_feat_atom(feat_coords)
            if atom is None:
                continue
            ideal_vector = feat.protein_vector
            coords = self._get_ideal_coordinates(ideal_vector, feat_coords, 3)
            if feat.type == 'donor_acceptor':
                d_score = self._get_near_score(coords, 'acceptor', tolerance=4)
                a_score = self._get_near_score(coords, 'donor', tolerance=4)
                self._update_score_dic(self.features_by_score, atom, feat, a_score, 'acceptor')
                for n in atom.neighbours:
                    if n.atomic_number == 1:
                        self._update_score_dic(self.features_by_score, n, feat, d_score, 'donor')

            elif atom_type_dict[feat.type] == 'acceptor':
                score = self._get_near_score(coords, atom_type_dict[feat.type], tolerance=3)
                for n in atom.neighbours:
                    if n.atomic_number == 1:
                        self._update_score_dic(self.features_by_score, n, feat, score, 'donor')
            elif atom_type_dict[feat.type] == 'donor':
                score = self._get_near_score(coords, atom_type_dict[feat.type], tolerance=3)
                self._update_score_dic(self.features_by_score, atom, feat, score, 'acceptor')

        return self.features_by_score

    def _score(self, values):
        '''Calculate geometric mean of scores'''
        return reduce(operator.__mul__, values, 1.0) ** (1. / len(values))

    def score_ligand_atoms(self, mol, dict=False, schematic=False, tolerance=2):
        '''
        Score ligand atoms in hotspot maps, taking interaction type into account.

        :param mol: a :class:`ccdc.Molecule` instance
        :param schematic: bool, catgorises score,for use with 2D depiction
        :param tolerance: int, the highest scoring grid point within +/- tolerance in each of the x,y and z directions will be assigned to the atom
        :return:

        '''
        if dict == True or schematic == True:
            scores = {}
        else:
            scores = []
        atom_labels = []
        for atom in mol.heavy_atoms:
            atom_labels.append(atom.label)
            atom_type = self._get_atom_type(atom)
            if atom_type == 'doneptor':
                d_score = self._get_near_score(atom.coordinates, 'donor', tolerance=tolerance)[0]
                a_score = self._get_near_score(atom.coordinates, 'acceptor', tolerance=tolerance)[0]
                score = max(d_score, a_score)
            else:
                score, coords = self._get_near_score(atom.coordinates, atom_type, tolerance=tolerance)
                # print(score, atom_type, atom.atomic_symbol)

            if schematic == True:
                if score >= 17:
                    score = 20
                elif score < 17 and score >= 14:
                    score = 10
                else:
                    score = 0
                scores.update({atom.label: score})

            elif dict == True:
                scores.update({str(atom.label): score})
            else:
                scores.append(score)
        # return atom_labels , scores
        return scores

    def schematic_map(self, ligand, title=False, output="diagram.png"):
        '''
        Display the distribution of scores as a heatmap on a 2D depiction of the molecule

        :param ligand: a :class:`ccdc.Molecule` object.
        :param title: str, Title placed at the top of the image
        :param output: str, Output file name
        :return:
        '''
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            from rdkit.Chem import AllChem
            from matplotlib.colors import LinearSegmentedColormap

        except ImportError:
            print("""rdkit is needed for this method""")
            exit()

        mol = MoleculeReader(ligand)[0]

        if ligand.split(".")[-1] == "mol2":
            with open(ligand, 'r') as lig:
                data = lig.read()
            mol_rdkit = Chem.MolFromMol2Block(data)
            AllChem.Compute2DCoords(mol_rdkit)

        elif ligand.split(".")[-1] == "sdf":
            suppl = Chem.SDMolSupplier(ligand)
            mol_rdkit = suppl[0]
            AllChem.Compute2DCoords(mol_rdkit)
        else:
            print("Method supports .mol2 files only!")
            raise ValueError

        scores = self.score_ligand_atoms(mol, schematic=True, tolerance=2)
        num_atoms = mol_rdkit.GetNumAtoms()
        a = 0.9 / (float(num_atoms))
        s = 0.005

        contribs = [float(scores[mol_rdkit.GetAtomWithIdx(i).GetProp('_TriposAtomName')]) for i in
                    range(mol_rdkit.GetNumAtoms())]

        fig = Draw.MolToMPL(mol_rdkit)

        cm = colourmap(scheme="inferno")
        test_cm = LinearSegmentedColormap.from_list(__file__, cm)
        try:
            x, y, z = Draw.calcAtomGaussians(mol_rdkit, a, step=s, weights=contribs)

            fig.axes[0].imshow(z, cmap=test_cm, interpolation='bilinear', origin='lower', alpha=0.9,
                               extent=(0, 1, 0, 1))
            fig.axes[0].contour(x, y, z, 5, colors='k', alpha=0.2)
        except ValueError:
            print("")

        if title:
            fig.text(1.25, 2.3, title, fontsize=20, horizontalalignment='center', verticalalignment='top',
                     color="white")

        fig.savefig(output, bbox_inches='tight')

    def score_ligand(self, mol, tolerance=2):
        '''
        Score ligand in hotspot maps, taking interaction type into account.

        :param mol: a :class:`ccdc.Molecule` object
        :return: float, Geometric mean of atomic scores
        '''
        scores = self.score_ligand_atoms(mol, dict=False, schematic=False, tolerance=tolerance)
        # print(scores)
        geo_mean = self._score(scores)
        return geo_mean

    def _get_percentage_rank(self, atom, atom_type):
        coordinates = atom.coordinates
        score, indices = self._get_near_score(coordinates, atom_type, 2)

        nx, ny, nz = self.super_grids[atom_type].nsteps

        g = self.super_grids[atom_type]
        scores = [g.value(i, j, k) for i in range(0, nx) for j in range(0, ny) for k in range(0, nz)
                  if g.value(i, j, k) > 0]

        filtered_scores = [x for x in scores if x < score]

        percentage_rank = 100 * float(len(filtered_scores)) / float(len(scores))
        return percentage_rank

    def _get_ligand_percentage_rank(self, mol):
        """
        Calculate the percentage rank for each atom of a ligand
        :param mol: a :class:`ccdc.Molecule` instance
        :return:
        """

        percentage_ranks = []
        for atom in mol.heavy_atoms:
            atom_type = self._get_atom_type(atom)
            if atom_type == 'doneptor':
                atom_type = 'donor'
            percentage_ranks.append(self._get_percentage_rank(atom, atom_type))
        return percentage_ranks

    def _get_grid_percentiles(self, g, percentiles):
        """
        percentile of non-zero grid values

        :param g:
        :param percentiles:
        :return:
        """

        nx, ny, nz = g.nsteps
        percentiles_dict = {}
        values = [g.value(i, j, k) for i in range(0, nx) for j in range(0, ny) for k in range(0, ny)
                  if g.value(i, j, k) > 1]
        if len(values) == 0:
            return percentiles_dict
        else:
            for strength, percentile in percentiles.items():
                cutoff = np.percentile(values, percentile) - 0.01
                percentiles_dict[strength] = cutoff

        return percentiles_dict

    def _single_grid(self):
        '''
        Takes all grids in the hotspot result object and applies a mask, such that for each grid point is only described
        by the highest scoring map. The grids are kept as separate object in order to track interaction type, however
        they could be summed to give a single grid
        :return:

        '''
        for probe, g in self.super_grids.items():
            masked_grids = {}
            for probe, g in self.super_grids.items():
                other_grids = [self.super_grids[p] for p in self.super_grids.keys() if p != probe]
                mg = g * (
                        (g > other_grids[0]) & (g > other_grids[1]))  # & (g >other_grids[2]) & (g >other_grids[3]))
                masked_grids[probe] = mg

        return masked_grids

    def _grid_values_in_range(self, grid, r_start, r_finish, score_cutoff):

        all_points = []
        range_points = []

        island = grid.islands(score_cutoff - 0.01)
        for g in island:
            nx, ny, nz = g.nsteps
            island_all_points = [g.value(i, j, k)
                                 for i in range(nx) for j in range(ny) for k in range(nz)
                                 if g.value(i, j, k) >= score_cutoff]

            island_range_points = [g.value(i, j, k)
                                   for i in range(nx) for j in range(ny) for k in range(nz)
                                   if g.value(i, j, k) >= r_start and g.value(i, j, k) < r_finish]

            all_points += island_all_points
            range_points += island_range_points

        avg_score = np.mean(all_points)
        return range_points, all_points, avg_score

    def _count_non_zero(self, cutoff):
        # print(cutoff)
        try:
            tmp_g = self._island_by_volume(self.sg, cutoff)

        except IndexError:
            return 9999
        try:
            nx, ny, nz = tmp_g.nsteps
        except AttributeError:
            return 9999
        points = [(i, j, k)
                  for i in range(nx) for j in range(ny) for k in range(nz)
                  if tmp_g.value(i, j, k) >= cutoff]

        return abs(self.num_gp - len(points))

    def _count_non_zero_alt(self, cutoff):
        print(cutoff)
        try:
            tmp_g = self._island_by_volume(self.sg, cutoff)

        except IndexError:
            return 999999
        try:
            nx, ny, nz = tmp_g.nsteps
        except AttributeError:
            return 0

        points = [(i, j, k)
                  for i in range(nx) for j in range(ny) for k in range(nz)
                  if tmp_g.value(i, j, k) >= cutoff]

        return self.num_gp - len(points)

    def _island_by_volume(self, sg, score_cutoff):
        '''

        :param sg:
        :param volume:
        :return:
        '''
        island = sg.islands(score_cutoff)
        island_by_total = {}
        for g in island:
            nx, ny, nz = g.nsteps
            island_points = [g.value(i, j, k)
                             for i in range(nx) for j in range(ny) for k in range(nz)
                             if g.value(i, j, k) >= score_cutoff]
            total = sum(island_points)
            island_by_total[total] = g

        try:
            top_score = sorted(island_by_total.keys(), reverse=True)[0]
            tmp_g = island_by_total[top_score]
        except IndexError:
            return None
        return tmp_g

    def _get_grid_by_volume(self, sg, volume):
        '''
        Takes a single grid input (produced by _single_grid()) and produces a mask that will provide a grid with the
        desired volume
        :param sg:
        :return:
        '''
        # Calculate number of grid points that corresponds to selected volume
        num_gp = int(float(volume) / 0.125)
        all_points = {}

        self.sg = sg
        self.num_gp = num_gp

        # Find the largest island at score_cutoff
        score_cutoff = optimize.fminbound(self._count_non_zero, 0, 30, xtol=0.025)

        if self._count_non_zero_alt(score_cutoff) < 0:
            score_cutoff -= 0.025
        if score_cutoff >= 29:
            score_cutoff = 1
        # print(score_cutoff)

        tmp_g = self._island_by_volume(self.sg, score_cutoff)

        if tmp_g is None:
            tmp_g = sg * 0
            return tmp_g
        c_sg, top_g = self._common_grid(sg, tmp_g, padding=0)

        nx, ny, nz = tmp_g.nsteps
        print(nx, ny, nz)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    coords = self._indices_to_point(i, j, k, tmp_g)
                    try:
                        all_points[c_sg.value_at_point(coords)].append(coords)
                    except KeyError:
                        all_points[c_sg.value_at_point(coords)] = [coords]

        scores = sorted(all_points.keys(), reverse=True)
        kept_points = []

        for score in scores:
            kept_points += all_points[score]
            new_cutoff = score
            if len(kept_points) >= num_gp:
                break
        print("Cutoff", new_cutoff)
        mask = (top_g > new_cutoff)
        return mask

    def _percent_in_range(self, grid, r_start, r_finish, score_cutoff):

        range_points, all_points, score = self._grid_values_in_range(grid, r_start, r_finish, score_cutoff)

        num_all = float(len(all_points))
        num_range = float(len(range_points))
        if num_all == 0:
            return 0

        return (num_range / num_all) * 100

    def _percent_by_type(self, grid_dict):

        all_g = grid_dict.values()
        sum_g = all_g[0].copy()

        for g in all_g[1:]:
            sum_g += g
        combined_range_points, combined_all_points, avg_score = self._grid_values_in_range(sum_g, 1, 10000, 1)

        percent_by_type = {}
        percent_by_type['average'] = avg_score
        percent_by_type['volume'] = float(len(combined_all_points)) * 0.125
        print(percent_by_type['volume'])
        for probe, g in grid_dict.items():

            range_points, all_points, score = self._grid_values_in_range(g, 1, 10000, 1)
            print(probe, len(combined_all_points), len(all_points))
            print(probe, len(combined_range_points), len(range_points))
            try:
                p = (float(len(all_points)) / float(len(combined_all_points))) * 100
            except ZeroDivisionError:
                p = 0
            percent_by_type[probe] = p

        return percent_by_type

    def _get_large_cavity(self, buriedness, drug_volume=300):
        """
        return cavities over a threshold size (volume approximately the volume of an average drug)

        :return:
        """

        large = []
        for island in buriedness.islands(threshold=5):
            x, y, z = island.nsteps
            size = len([island.value(nx, ny, nz) for nx in range(x) for ny in range(y) for nz in range(z)
                        if island.value(nx, ny, nz) > 5]) * 0.125
            if size > drug_volume:
                large.append(island)
        return large

    def drug_sized(self, buriedness):
        """
        detemines whether a hotspot resides within a region that can accomodate a drug sized molecule

        :return:
        """

        bottom_left, top_right = self.super_grids["apolar"].bounding_box
        centroid = molecule.Coordinates(((bottom_left.x + top_right.x) / 2),
                                        ((bottom_left.y + top_right.y) / 2),
                                        ((bottom_left.z + top_right.z) / 2)
                                        )
        drug_cavity = self._get_large_cavity(buriedness=buriedness)

        for d in drug_cavity:
            if self.contains_point(centroid, d):
                return True
            else:
                return False

    def hotspot_score(self, total=True):
        """
        caution! only use for extracted hotspots!
        :return:
        """

        if total:
            t = [self.grid_score(g, threshold=0) for probe, g in self.super_grids.items()
                 if self.grid_score(g, threshold=0) > 0]
            return np.mean(t)

        else:
            return self.grid_score(self.super_grids["apolar"], threshold=0)

    def get_bfactors(self, pdb_file):
        """
        Assign B factors
        """

        lines = self._get_lines_from_file(pdb_file)

        return {str(l[21] + ":" + l[17:20] + l[23:26] + "_" + l[12:16]).replace(" ", ""): float(l[60:66])
                for i, l in enumerate(lines) if str(l[0:4]) == "ATOM"}

    def hotspot_bfactor(self, pdb_file, bfactor):
        """


        :param prot:
        :param bfactor:
        :return:
        """
        prot = Protein.from_file(pdb_file)
        prot.remove_all_metals()
        prot.remove_all_waters()
        for lig in prot.ligands:
            prot.remove_ligand(lig.identifier)

        mini, maxi = self.super_grids["apolar"].bounding_box
        centroid = molecule.Coordinates((mini.x + maxi.x) / 2, (mini.y + maxi.y) / 2, (mini.z + maxi.z) / 2)
        bs = prot.BindingSiteFromPoint(protein=prot, origin=centroid, distance=12)
        residues = [str(r.identifier) for r in bs.residues]
        protein = [str(p.identifier) for p in prot.residues]
        deletes = list(set(protein) - set(residues))
        for delete in deletes:
            prot.remove_residue(delete)
        b = [bfactor[str(r.identifier + "_" + a.label)] for r in prot.residues for a in r.atoms]
        self.prot = prot
        print("b", b)
        if len(b) == 0:
            return None, None
        else:
            return np.percentile(b, 50), np.mean(b)

    def output_data(self, hr, buriedness, pdb_file):
        """

        :param hr:
        :return:

        """

        identifier = []
        donor = []
        acceptor = []
        negative = []
        positive = []
        drug_sized = []
        apolar_score = []
        hotspot_score = []
        bfactor_mean = []
        bfactor_median = []
        hotspot_bfactor_mean = []
        hotspot_bfactor_median = []

        bfactor = self.get_bfactors(pdb_file)
        for i, h in enumerate(hr):
            composition = {probe: len(g.islands(threshold=14)) for probe, g in h.super_grids.items()}
            identifier.append(i)
            donor.append(composition["donor"])
            acceptor.append(composition["acceptor"])
            negative.append(composition["negative"])
            positive.append(composition["positive"])
            drug_sized.append(h.drug_sized(buriedness))
            apolar_score.append(h.hotspot_score(total=False))
            hotspot_score.append(h.hotspot_score(total=True))
            bfactor_mean.append(np.mean(bfactor.values()))
            bfactor_median.append(np.percentile(bfactor.values(), 50))
            mean, median = h.hotspot_bfactor(pdb_file, bfactor)
            if pdb_file is not None:
                hotspot_bfactor_mean.append(mean)
                hotspot_bfactor_median.append(median)

        return pd.DataFrame({"identifier": identifier, "donor": donor, "acceptor": acceptor, "negative": negative,
                             "positive": positive, "drug_sized_cavity": drug_sized, "apolar_av_score": apolar_score,
                             "average_score": hotspot_score, "bf_mean": bfactor_mean, "bf_median": bfactor_median,
                             "hbf_mean": hotspot_bfactor_mean, "hbf_median": hotspot_bfactor_median
                             })

    def extract_hotspots(self, out_dir=None, pharmacophore=True, **kwargs):
        """
        For a given output volume, hotspots are identified by the peaks in apolar propensity.

        :param out_dir:
        :param sigma: float, target volume to be selected in Angstroms^3
        :param cutoff: int, threshold value to contour islands
        :param volume: int, volume of the desired output hotspots.
        :return:
        """
        fragments = kwargs.get("fragment")
        lead = kwargs.get("lead")
        bfactors = kwargs.get("bfactors")

        build = HotspotBuilder(self.super_grids, self.sampled_probes, out_dir, self.prot, kwargs)
        build.write(mode="locations")          # for testing
        build.write(mode="features")           # for testing
        build.construct()
        build.format_data()

        # create hotspot result
        hrs = []
        for i, hotspot in enumerate(build.locations):
            print(hotspot.extracted_super_grids)
            hotspot.identifier = str(i)
            hr = HotspotResults(grid_dict=hotspot.extracted_super_grids,
                                protein=self.prot,
                                fname=self.fname,
                                sampled_probes=hotspot.extracted_probes,
                                buriedness=self.buriedness,
                                out_dir=self.out_dir)
            hrs.append(hr)

            # create pharmacophore
            if pharmacophore:
                hr.pharmacophore = hr.get_pharmacophore_model(identifier=hotspot.identifier)
                hr.pharmacophore.identifier = hotspot.identifier
                #hr.pharmacophore.fname = pharmacophore
                hrs.append(hr)

        build.pymol_out = extracted_hotspot_template(len(hrs), fragments, lead)
        for hr in hrs:
            if hr.pharmacophore is not None:
                build.pymol_out += hr.pharmacophore.get_pymol_pharmacophore()

        with open(join(out_dir, "hotspot_boundaries","extracted_hotspots.py"), 'w') as w:
            w.write(build.pymol_out)

        with MoleculeWriter(join(out_dir, "hotspot_boundaries", "protein.pdb")) as w:
            w.write(self.prot)

        #self.output_extracted_hotspots(hrs, out_dir, fragments, lead, charged=False)
        #
        # if bfactors is not None:
        #     df = self.output_data(hrs, self.buriedness, bfactors)
        #     df.to_csv(join(out_dir, "data.csv"))

        return hrs

    def best_continuous_volume(self, volume=500, pocket_mask=False):
        '''
        Selects the highest scoring continuous region of propensity across all map types, matching the input volume.
        Selecting a molecular volume corresponding to a typical drug-like molecule (350-550 Angstroms ^3) provides a
        prediction of the most tractable binding site on the protein.

        :param volume: float, target volume to be selected in Angstroms^3
        :param pocket_mask: a :class:`ccdc.interaction.Grid` instance with a value of 1 at any grid points to be ignored, otherwise 0.
        :return:
        '''

        processed_grids = {}
        # print("start descriptors")

        masked_grids = self._single_grid()
        all_g = masked_grids.values()
        sum_g = all_g[0].copy()
        for g in all_g[1:]:
            sum_g += g

        if pocket_mask:
            com_sum_grid, com_pocket_mask = self._common_grid(sum_g, pocket_mask, padding=0)
            sum_g = com_sum_grid * (com_pocket_mask < 1)

        second_mask = self._get_grid_by_volume(sum_g, volume)

        for probe, g in self.super_grids.items():
            print(probe, )
            mg = masked_grids[probe]
            c_mg, c_second_mask = self._common_grid(mg, second_mask, padding=0)
            out_g = (c_mg * c_second_mask)

            del (c_mg)
            del (c_second_mask)

            # out_g.write('processed_{}.grd'.format(probe))
            processed_grids[probe] = out_g

        bcv_hr = HotspotResults(processed_grids, self.prot, self.fname, None, None,
                                out_dir='{}_bcv'.format(self.out_dir))

        remaining = {}
        for probe, g in self.super_grids.items():
            diff_g = g - bcv_hr.super_grids[probe]
            remaining.update({probe: diff_g})

        remaining_hr = HotspotResults(remaining, self.prot, self.fname, None, None,
                                      out_dir='{}_remaining'.format(self.out_dir))

        return bcv_hr, remaining_hr

    def extract_pocket(self, whole_residues=False):
        '''
        Create a :class:`ccdc.Protein` containing atoms or residues that have a score

        :param whole_residues: bool, whether to include all residue atoms if only a subset have a score > 0
        :return: a :class:`ccdc.Protein` instance
        '''
        prot_scores = self.score_protein()
        pocket = self.prot.copy()
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
    
    def _smooth_crude(self):
        for probe, g in self.super_grids.items():
            self.super_grids[probe] = g.max_value_of_neighbours()

    def zip_results(self, archive_name='out', delete_directory = True):
        """
        Zips the output directory created for this :class:`fragment_hotspot_maps.HotspotResults` instance, and
        removes the directory by default. The zipped file can be loaded directly into a new
        :class:`fragment_hotspot_maps.HotspotResults` instance using the 
        :func:`~fragment_hotspot_maps.Hotspots.from_zip_dir` function

        :param archive_name: str, file path
        :param delete_directory: bool, remove the out directory once it has been zipped

        :return: None
        """
        print(getcwd())
        if not exists(self.out_dir):
            mkdir(self.out_dir)

        for probe, g in self.super_grids.items():
            g.write('{0}/{1}.grd'.format(self.out_dir, probe))

        with MoleculeWriter(join(self.out_dir, "protein.pdb")) as writer:
            writer.write(self.prot)

        self.archive_name = archive_name
        shutil.make_archive(self.archive_name, 'zip', self.out_dir)
        self.archive_loc = dirname("{}.zip".format(self.archive_name))
        if delete_directory:
            shutil.rmtree(self.out_dir)

        

    def output_pymol_file(self):
        """
        Output a python script to be run from within pymol to visualise output
        :return:
        """

        percentiles = {'low': 30}
        cutoff_by_probe = {}
        for probe, g in self.super_grids.items():
            cutoff_dict = self._get_grid_percentiles(g, percentiles)
            cutoff_by_probe[probe] = cutoff_dict.values()

        if self.archive_name is not None:

            with open(path.join(self.archive_loc,'pymol_results_file_{}.py'.format(self.archive_name)), 'w') as pymol_file:
                pymol_out = pymol_zip_template(self.archive_name, self.super_grids.keys(), cutoff_by_probe)
                pymol_file.write(pymol_out)

        else:

          
            if self.fname is None:
                with MoleculeWriter('{}/protein.pdb'.format(self.out_dir)) as w:
                    w.write(self.prot)
                    prot_file = ('{}/protein.pdb'.format(self.out_dir))
            else:
                prot_file = self.fname
    
    
            for probe, g in self.super_grids.items():
                g.write('{0}/{1}.grd'.format(self.out_dir, probe))
    
            with open('{}/pymol_results_file.py'.format(self.out_dir), 'w') as pymol_file:
                pymol_out = pymol_template(prot_file, self.out_dir, self.super_grids.keys(), cutoff_by_probe)
                pymol_file.write(pymol_out)
    
            with MoleculeWriter(join(self.out_dir, "protein.pdb")) as writer:
                writer.write(self.prot)


class Hotspots(_HotspotsHelper):
    """
    A class for running Fragment Hotspot Map calculations
    """

    class _Sampler(_HotspotsHelper):
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
            self.probe_grids = [_SampleGrid(g.name, self._copy_and_clear(g.grid), g.atom_predicate) for g in self.grids]

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
                    d = self._get_distance(c, a.coordinates)
                    atom_by_distance[d] = a
            else:
                for a in molecule.atoms:
                    d = self._get_distance(c, a.coordinates)
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

                maxima = [self._indices_to_point(i, j, k, g)
                          for i in range(nx) for j in range(ny) for k in range(nz)
                          if g.value(i, j, k) >= translation_threshold]

                translate_probe = translate_probe + maxima
            # print(priority_atom_type, len(translate_probe))
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
                    i, j, k = self._point_to_indices(pg.add_coordinates(a, trans), pg.grid)
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

            print("Q", len(quaternions), "P", probe)

            for g in self.grids:
                g.set_molecule(molecule, True)

            for g in self.probe_grids:
                g.set_molecule(molecule, self.settings.polar_contributions)

            for q in quaternions:
                molecule.apply_quaternion(q)
                priority_atom_coordinates = priority_atom.coordinates
                active_coordinates_dic = self.get_active_coordinates()

                for priority_atom_point in translate_points:

                    translation = [priority_atom_point[i] - priority_atom_coordinates[i]
                                   for i in xrange(0, len(priority_atom_coordinates))]

                    score = self.sample_pose(translation, active_coordinates_dic, probe)

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
                    self.update_out_grids(score, active_coordinates_dic, translation)

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

        self.superstar_grids = None
        self.weighted_grids = None
        self.sampled_probes = {}

        self.prot = None
        self.fname = None
        self.probe_size = None
        self.charged_probes = None

        self.out_dir = None
        self.wrk_dir = None

        if settings is None:
            self.sampler_settings = self._Sampler.Settings()
        else:
            self.sampler_settings = settings

    def _get_cavity_centroid(self, cav):
        """
        Returns the centroid of a cavity object

        :param cav:
        :return:
        """

        x_coords = []
        y_coords = []
        z_coords = []

        for feat in cav.features:
            feature_coords = feat.coordinates
            x_coords.append(feature_coords[0])
            y_coords.append(feature_coords[1])
            z_coords.append(feature_coords[2])

        x_avg = round(np.mean(x_coords))
        y_avg = round(np.mean(y_coords))
        z_avg = round(np.mean(z_coords))

        return x_avg, y_avg, z_avg

    @staticmethod
    def _superstar_job(args):
        """
        creates a RunSuperstar class and returns a SuperstarResult object
        :param args: tuple, probe name and centroid
        :return: a :class: `__main__().SuperstarResult` object
        """

        s = _RunSuperstar()
        n, ss_probe, centroid, prot, out_dir, wrk_dir = args
        s.settings.jobname = "{}.ins".format(n)
        s.settings.probename = ss_probe
        s.settings.moleculefile = "protein.pdb"
        s.settings.cavity_origin = centroid
        s.settings.working_directory = wrk_dir
        result = s.run_superstar(prot, out_dir)
        return result

    def _run_ss(self, centroid=None):
        """
        initiates a SuperStar run for a given protein and probe

        :param prot: a :class:`ccdc.protein.Protein` instance
        :param out_dir: str, output directory
        :param centroid: tup, coordinates of cavity origin
        :param charged_probes: bool, if True 'positive' and 'negative' probes will be used
        :return: a :class:`SuperstarResult` instance
        """

        if self.charged_probes:
            self.probe_dict = dict(
                apolar='AROMATIC CH CARBON',
                donor='UNCHARGED NH NITROGEN',
                acceptor='CARBONYL OXYGEN',
                positive='CHARGED NH NITROGEN',
                negative='CARBOXYLATE OXYGEN'
            )
        else:
            self.probe_dict = dict(
                apolar='AROMATIC CH CARBON',
                donor='UNCHARGED NH NITROGEN',
                acceptor='CARBONYL OXYGEN')

        wrk_dir = _test_output_dir()
        args = [(k, self.probe_dict[k], centroid, self.prot, self.out_dir, wrk_dir) for k in self.probe_dict.keys()]
        ex = futures.ThreadPoolExecutor(max_workers=5)
        results = ex.map(self._superstar_job, args)

        return list(results)

    def _get_weighted_maps(self):
        """
        weight superstar output by burriedness

        :return: a list of :class: `WeightedResult` instances
        """

        results = []
        for s in self.superstar_grids:
            if self.ghecom_executable:
                self.buriedness = self.ghecom.grid
            else:
                self.buriedness = s.ligsite

            s.grid, self.buriedness = self._common_grid(s.grid, self.buriedness, 1)

            weighted_grid = s.grid * self.buriedness
            results.append(_WeightedResult(s.identifier, weighted_grid))
        return results

    def _get_out_maps(self, probe, grid_dict):
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

        probe_path = pkg_resources.resource_filename('fragment_hotspot_maps', 'probes/')

        if self.charged_probes:
            if probe == "negative" or probe == "positive":
                mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, "test")))[0]
            else:
                mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, self.probe_size)))[0]
        else:
            mol = MoleculeReader(join(probe_path, "rotate-{}_{}_flat.mol2".format(probe, self.probe_size)))[0]

        probes = self.sampler.sample(mol, probe=probe)
        for pg in self.sampler.probe_grids:
            if pg.name.lower() == probe:
                try:
                    self.out_grids[pg.name].append(pg.grid)
                except KeyError:
                    self.out_grids[pg.name] = [pg.grid]
        return probes

    def _calc_hotspots(self):
        """
        Function for overall organisation of hotspot calculation

        :return:
        """

        print("Start SS")

        probe_types = ['apolar', 'donor', 'acceptor']
        if self.charged_probes:
            probe_types += ['negative', 'positive']

        self.superstar_grids = self._run_ss(self.centroid)

        print("SS complete")

        if self.ghecom_executable:
            out_grid = self._copy_and_clear(self.superstar_grids[0].ligsite)
            r = _RunGhecom()
            r.settings.prot = self.prot
            r.settings.out_grid = out_grid
            r.settings.ghecom_executable = self.ghecom_executable
            self.ghecom = r.run_ghecom()
            chdir(self.working_dir)


        self.weighted_grids = self._get_weighted_maps()
        self.buriedness.write(join(self.out_dir, "buriedness.grd"))

        grid_dict = {w.identifier: w.grid for w in self.weighted_grids}

        for probe in probe_types:
            top_probes = self._get_out_maps(probe, grid_dict)
            self.sampled_probes.update({probe: top_probes})

    def _from_grid_dic(self, super_grids, prot, fname=None, sampled_probes=None, buriedness=None):
        """
        Create a Hotspots_reults object from a dictionary of previously calculated grid objects

        :param super_grids: dict, {'probe_name':grid} where the probe names are 'apolar', 'donor' and 'acceptor'
        :param prot: a :class:`ccdc.protein.Protein` instance
        :param fname: str, file path
        :return: a :class:`fragment_hotspot_maps.HotspotResults` instance
        """

        self.fname = fname
        self.super_grids = super_grids
        self.prot = prot
        self.sampled_probes = sampled_probes
        self.buriedness = buriedness
        return HotspotResults(self.super_grids, self.prot, self.fname, self.sampled_probes, self.buriedness,
                              out_dir=self.out_dir)

    def from_zip_dir(self, zip_dir, probes = ('apolar','donor','acceptor'), output_directory= 'out', temporary = True):
        """
        Create a Hotspots_reults object from zipped output directory with default file names. Extracts files to a
        temp dir, and removes files once a :class:`fragment_hotspot_maps.HotspotResults` instance has been created

        :param super_grids: dict, {'probe_name':grid} where the probe names are 'apolar', 'donor' and 'acceptor'
        :param prot: a :class:`ccdc.protein.Protein` instance
        :param fname: str, file path
        :return: a :class:`fragment_hotspot_maps.HotspotResults` instance
        """

        self.fname = 'protein.pdb'
        self.out_dir = output_directory

        if temporary:
            dirpath = tempfile.mkdtemp()
        else:
            dirpath = self.out_dir
            if not exists(self.out_dir):
                mkdir(self.out_dir)


        with zipfile.ZipFile(zip_dir) as hs_zip:
            hs_zip.extractall(dirpath)

        for p in probes:
            print(path.join(dirpath,'{}.grd'.format(p)))

        self.super_grids = {p: Grid.from_file(path.join(dirpath,'{}.grd'.format(p))) for p in probes}
        self.prot = Protein.from_file(path.join(dirpath,'protein.pdb'))
        self.sampled_probes = {}
        try:
            self.buriedness = Grid.from_file(path.join(dirpath,'buriedness.grd'))
        except:
            self.buriedness=None

        if temporary:
            shutil.rmtree(dirpath)
        return HotspotResults(self.super_grids, self.prot, self.fname, self.sampled_probes, self.buriedness,
                              out_dir=self.out_dir)

    def from_protein(self, prot, charged_probes=False, binding_site_origin=None, probe_size=7,
                     ghecom_executable=None, output_directory = 'out'):
        """
        Calculate Fragment Hotspot Maps from a ccdc.protein.Protein object

        :param prot: a :class:`ccdc.protein.Protein` instance
        :param binding_site_origin: a tuple of three floats, giving a coordinate within the binding site
        :param probe_size: int, size of probe in number of heavy atoms (3-8 atoms)
        :param ghecom_executable: str, path to ghecom executeable, if None ligsite used
        :param charged_probes: bool, if True include positive and negative probes
        :return: a :class:`fragment_hotspot_maps.HotspotResults` instance
        """

        self.prot = prot
        self.fname = 'protein.pdb'
        self.probe_size = probe_size
        self.charged_probes = charged_probes
        self.out_dir = output_directory
        if not exists(self.out_dir):
            mkdir(self.out_dir)
        self.out_dir = path.abspath(self.out_dir)
        self.working_dir = getcwd()

        self.ghecom_executable = ghecom_executable
        self.centroid = binding_site_origin

        with MoleculeWriter(join(self.out_dir, 'protein.pdb')) as w:
            w.write(self.prot)

        self._calc_hotspots()

        for probe in self.out_grids.keys():
            probe = probe.lower()
            sg = self.out_grids[probe][0]
            self.super_grids[probe] = sg
        return HotspotResults(self.super_grids, self.prot, self.fname, self.sampled_probes, self.buriedness,
                              out_dir=self.out_dir)


def main():
    prot_file = "Z:/fragment_hotspot_maps/protein.pdb"
    ghecom_exe = None
    prot = Protein.from_file(prot_file)

    h = Hotspots()

    result = h.from_protein(prot=prot, charged_probes=False, fname=prot_file, ghecom_executable=ghecom_exe)
    result.output_pymol_file(out_dir=None)


if __name__ == "__main__":
    main()
