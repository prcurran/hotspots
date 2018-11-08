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

The :mod:`fragment_hotspot_maps.utilities` module contains classes to for
general functionality.

The main classes of the :mod:`fragment_hotspot_maps.extraction` module are:
    -Utilities
"""
import collections
import math
from os.path import exists, join
from os import mkdir
import tempfile

import numpy as np

from ccdc.io import MoleculeWriter
from ccdc.cavity import Cavity

Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])

class Utilities(object):
    """
    class providing utility functions
    """

    @staticmethod
    def get_distance(coords1, coords2):
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
    def get_out_dir(path):
        """"""
        if not exists(path):
            mkdir(path)
        return path

    @staticmethod
    def get_lines_from_file(fname):
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

    @staticmethod
    def cavity_centroid(obj):
        """
        Returns the centroid of a cavity object
        :param cav:
        :return:
        """
        x_coords = []
        y_coords = []
        z_coords = []
        if isinstance(obj, Cavity):
            features = (f.coordinates for f in cav.features)

        else:
            features = obj.surface_points

            for feat in features:
                x_coords.append(feat[0])
                y_coords.append(feat[1])
                z_coords.append(feat[2])

        x_avg = round(np.mean(x_coords))
        y_avg = round(np.mean(y_coords))
        z_avg = round(np.mean(z_coords))

        return Coordinates(x=x_avg, y=y_avg, z=z_avg)

    @staticmethod
    def cavity_from_protein(prot):
        """
        generates a cavities from prot.
        :param prot:
        :return:
        """
        # TODO: proper solution has been requested

        tfile = join(tempfile.mkdtemp(), "protein.pdb")
        with MoleculeWriter(tfile) as writer:
            writer.write(prot)

        return Cavity.from_pdb_file(tfile)

class Figures(object):
    """
    Class to handle the generation of hotspot related figures

    TO DO: is there a better place for this to live?
    """
    def histogram(self, hr, plot):
        """ creates a histrogram from grid values"""
        data = {}
        for p, g in hr.super_grids:
            array = g.get_array()
            masked_array = ma.masked_less_equal(array, 1)
            grid_values = masked_array.compressed()
            data.update({p, grid_values})

        if plot:
            plt = self._plot_histogram(data)
            return data, plt
        else:
            return data

    @staticmethod
    def _2D_diagram(hr, ligand, fpath, title):
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

        scores = hr.score_ligand_atoms(mol, schematic=True, tolerance=2)
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


    def _plot_histogram(self, data):
            """
            initialise the matplotlib figure to output histograms

            :param data:
            :return:
            """
            plt.figure(1)
            for n, key in enumerate(data.keys()):
                plt.subplot(3, 1, (n + 1))
                hist, bin_edges = np.histogram(data[key], bins=range(0, 40), normed=True)
                self._histogram_settings(bin_edges)
                if n == 0:
                    plt.title(title)
                if n < 2:
                    plt.xticks([])
                if n == 2:
                    plt.xlabel("Fragment hotspot score")
                if n == 1:
                    plt.ylabel("Frequency")
                plt.bar(bin_edges[:-1], hist, width=1, color=colour_dict[key])
            return plt

    @staticmethod
    def _histogram_settings(bin_edges):
        """
        display settings
        :param bin_edges:
        :param title:
        :return:
        """
        colour_dict = {"acceptor": "r", "donor": "b", "apolar": "y"}
        plt.xlim(min(bin_edges), max(bin_edges))
        plt.ylim(0, 0.35)
        plt.yticks([])

# def _histogram_info(self, data, key, n):
#     """
#     data for the matplotlib figure
#
#     :param data:
#     :param key:
#     :param n:
#     :return:
#     """
#
#     colour_dict = {"acceptor": "r", "donor": "b", "apolar": "y"}
#     hist, bin_edges = np.histogram(data[key], bins=range(0, 40), normed=True)
#     plt.bar(bin_edges[:-1], hist, width=1, color=colour_dict[key])
#     plt.xlim(min(bin_edges), max(bin_edges))
#     plt.ylim(0, 0.35)
#     plt.yticks([])
#     if n == 0:
#         plt.title("Fragment-hotspot Maps")
#     if n < 2:
#         plt.xticks([])
#     if n == 2:
#         plt.xlabel("Fragment hotspot score")
#     if n == 1:
#         plt.ylabel("Frequency")
#
#
# def _generate_histogram(self, data):
#     """
#     initialise the matplotlib figure to output histograms
#
#     :param data:
#     :return:
#     """
#
#     plt.figure(1)
#     for n, key in enumerate(data.keys()):
#         plt.subplot(3, 1, (n + 1))
#         self._histogram_info(data, key, n)
#     plt.savefig("hotspot_histogram.png")
#     plt.close()
#
#
# def hotspot_histograms(self, ligand=None):
#     """
#     Outputs histograms of hotspot score.
#
#     :param ligand:
#     :return: None
#     """
#
#     self.data = {}
#     if ligand:
#         self.bs_centroid = ligand.centre_of_geometry()
#         for g in self.super_grids.keys():
#             grd = self.super_grids[g]
#             fragment_centroid = self._point_to_indices(self.bs_centroid, grd)
#             n = 8
#             rx = range(fragment_centroid[0] - n, fragment_centroid[0] + n)
#             ry = range(fragment_centroid[1] - n, fragment_centroid[1] + n)
#             rz = range(fragment_centroid[2] - n, fragment_centroid[2] + n)
#
#             self.data[g] = np.array(
#                 [grd.value(i, j, k) for i in rx for j in ry for k in rz if grd.value(i, j, k) != 0])
#
#     else:
#         for g in self.super_grids.keys():
#             grd = self.super_grids[g]
#             nx, ny, nz = grd.nsteps
#             self.data[g] = np.array(
#                 [grd.value(i, j, k) for i in xrange(0, nx) for j in xrange(0, ny) for k in xrange(0, nz) if
#                  grd.value(i, j, k) != 0])
#
#     self._generate_histogram(self.data)
