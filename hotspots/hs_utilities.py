"""
The :mod:`hotspots.utilities` module contains classes to for
general functionality.

The main classes of the :mod:`hotspots.extraction` module are:
    - :class:`hotspots.hs_utilities.Helper`
    - :class:`hotspots.hs_utilities.Figures`

"""
from __future__ import division

import collections
import math
import tempfile
from os import mkdir
from os.path import exists, join, abspath

# import matplotlib as mpl
# mpl.use('TkAgg')
import matplotlib.pyplot as plt

import numpy as np

from ccdc.cavity import Cavity
from ccdc.io import MoleculeWriter
from ccdc.molecule import Molecule, Atom


Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


def _generate_usr_moment(fcoords_list):
    """
    PRIVATE EXPERIMENTAL FEATURE

    Modification of Adrian Schreyer's code https://bitbucket.org/aschreyer/usrcat

    More information on the USR algorithm:
        - https://doi.org/10.1098/rspa.2007.1823

    :param grds:

    :return:
    """
    from usrcat.geometry import usr_moments, usr_moments_with_existing

    all_coords = np.concatenate([c for c in fcoords_list if len(c) != 0])
    #np.vstack(fcoords_list)
    (ctd, cst, fct, ftf), om = usr_moments(all_coords)

    for fcoords in fcoords_list:

        # initial zeroed out USRCAT feature moments
        fm = np.zeros(12)

        # only attempt to generate moments if there are enough atoms available!
        if len(fcoords):
            fm = usr_moments_with_existing(fcoords, ctd, cst, fct, ftf)

        # append feature moments to the existing ones
        om = np.append(om, fm)

    return np.array([om])



class Helper(object):
    """
    A class to handle miscellaneous functionality
    """
    @staticmethod
    def get_distance(coords1, coords2):
        """
        given two coordinates, calculates the distance

        :param tup coords1: float(x), float(y), float(z), coordinates of point 1
        :param tup coords2: float(x), float(y), float(z), coordinates of point 2
        :return: float, distance
        """
        xd = coords1[0] - coords2[0]
        yd = coords1[1] - coords2[1]
        zd = coords1[2] - coords2[2]
        d = math.sqrt(xd ** 2 + yd ** 2 + zd ** 2)
        return d

    @staticmethod
    def get_out_dir(path):
        """
        checks if directory exists, if not, it create the directory

        :param str path: path to directory
        :return str: path to output directory
        """
        if not exists(abspath(path)):
            mkdir(abspath(path))
        return abspath(path)

    @staticmethod
    def get_lines_from_file(fname):
        """
        gets lines from text file, used in Ghecom calculation

        :return: list, list of str
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
        returns the centre of a cavity

        :param obj: can be a `ccdc.cavity.Cavity` or
        :return: Coordinate
        """
        if isinstance(obj, Cavity):
            features = [f.coordinates for f in obj.features]

        else:
            features = obj.surface_points

        x_avg = round(np.mean([feat[0] for feat in features if isinstance(feat[0], float)]))
        y_avg = round(np.mean([feat[1] for feat in features if isinstance(feat[1], float)]))
        z_avg = round(np.mean([feat[2] for feat in features if isinstance(feat[2], float)]))

        return Coordinates(x=x_avg, y=y_avg, z=z_avg)

    @staticmethod
    def cavity_from_protein(prot):
        """
        currently the Protein API doesn't support the generation of cavities directly from the Protein instance
        this method handles the tedious writing / reading

        :param `ccdc.protein.Protein` prot: protein
        :return: `ccdc.cavity.Cavity`
        """
        tfile = join(tempfile.mkdtemp(), "protein.pdb")
        with MoleculeWriter(tfile) as writer:
            writer.write(prot)

        return Cavity.from_pdb_file(tfile)

    @staticmethod
    def get_label(input, threshold=None):
        """
        creates a value labels from an input grid dictionary

        :param dic input: key = "probe identifier" and value = `ccdc.utilities.Grid`
        :return `ccdc.molecule.Molecule`: pseduomolecule which contains score labels
        """
        min_size_dict = {"apolar": 40,
                         "donor": 15,
                         "acceptor": 15,
                         "negative": 15,
                         "positive": 15}

        atom_dic = {"apolar": 'C',
                    "donor": 'N',
                    "acceptor": 'O',
                    "negative": 'S',
                    "positive": 'H'}

        try:
            interaction_types = [atom_dic[feat.feature_type] for feat in input._features]
            coordinates = [feat.feature_coordinates for feat in input._features]
            scores = [feat.score_value for feat in input._features]

        except AttributeError:

            print(threshold)
            try:
                if threshold is None:
                    pass
                else:
                    interaction_types = []
                    coordinates = []
                    scores = []
                    for p, g in input.items():
                        for island in g.islands(threshold=threshold):
                            if island.count_grid() > min_size_dict[p]:
                                interaction_types.append(atom_dic[p])
                                coordinates.append(island.centroid())
                                scores.append(max(island.grid_values(threshold=threshold)))
            except AttributeError:
                print("object not supported")

        mol = Molecule(identifier="pharmacophore_model")

        pseudo_atms = [Atom(atomic_symbol=interaction_types[i],
                            atomic_number=14,
                            coordinates=coordinates[i],
                            label=str(scores[i]))
                       for i in range(len(interaction_types))]

        for a in pseudo_atms:
            mol.add_atom(a)
        return mol

    @staticmethod
    def get_atom_type(atom):
        """
        return the atom classification

        :param atom:
        :return:
        """
        if atom.is_donor and \
                atom.atomic_symbol == 'N' and \
                len([a for a in atom.neighbours if a.atomic_symbol == 'H']) >= 2:
            return "donor"
        elif atom.is_donor and atom.is_acceptor:
            return "doneptor"
        elif atom.is_acceptor:
            return "acceptor"
        elif atom.is_donor:
            return "donor"
        else:
            return "apolar"


class Figures(object):
    """
    Class to handle the generation of hotspot related figures

    TO DO: is there a better place for this to live?
    """
    @staticmethod
    def histogram(hr):
        """
        creates a histogram from the hotspot scores

        :param `hotspots.result.Results` hr: a Fragment Hotspot Map result
        :return: data, plot
        """
        data = {}
        for p, g in hr.super_grids.items():
            array = g.get_array()
            masked_array = np.ma.masked_less_equal(array, 1)
            grid_values = masked_array.compressed()
            data.update({p: grid_values})

        plt = Figures._plot_histogram(data)
        return data, plt


    # @staticmethod
    # def _2D_diagram(hr, ligand, fpath, title):
    #     '''
    #     Display the distribution of scores as a heatmap on a 2D depiction of the molecule
    #
    #     :param ligand: a :class:`ccdc.Molecule` object.
    #     :param title: str, Title placed at the top of the image
    #     :param output: str, Output file name
    #     :return:
    #     '''
    # try:
    #     from rdkit import Chem
    #     from rdkit.Chem import Draw
    #     from rdkit.Chem import AllChem
    #     from matplotlib.colors import LinearSegmentedColormap
    # except ImportError:
    #     print("""rdkit is needed for this method""")
    #     exit()
    #
    #     mol = MoleculeReader(ligand)[0]
    #
    #     if ligand.split(".")[-1] == "mol2":
    #         with open(ligand, 'r') as lig:
    #             data = lig.read()
    #         mol_rdkit = Chem.MolFromMol2Block(data)
    #         AllChem.Compute2DCoords(mol_rdkit)
    #
    #     elif ligand.split(".")[-1] == "sdf":
    #         suppl = Chem.SDMolSupplier(ligand)
    #         mol_rdkit = suppl[0]
    #         AllChem.Compute2DCoords(mol_rdkit)
    #     else:
    #         print("Method supports .mol2 files only!")
    #         raise ValueError
    #
    #     scores = hr.score_ligand_atoms(mol, schematic=True, tolerance=2)
    #     num_atoms = mol_rdkit.GetNumAtoms()
    #     a = 0.9 / (float(num_atoms))
    #     s = 0.005
    #
    #     contribs = [float(scores[mol_rdkit.GetAtomWithIdx(i).GetProp('_TriposAtomName')]) for i in
    #                 range(mol_rdkit.GetNumAtoms())]
    #
    #     fig = Draw.MolToMPL(mol_rdkit)
    #
    #     cm = colourmap(scheme="inferno")
    #     test_cm = LinearSegmentedColormap.from_list(__file__, cm)
    #     try:
    #         x, y, z = Draw.calcAtomGaussians(mol_rdkit, a, step=s, weights=contribs)
    #
    #         fig.axes[0].imshow(z, cmap=test_cm, interpolation='bilinear', origin='lower', alpha=0.9,
    #                            extent=(0, 1, 0, 1))
    #         fig.axes[0].contour(x, y, z, 5, colors='k', alpha=0.2)
    #     except ValueError:
    #         print("")
    #
    #     if title:
    #         fig.text(1.25, 2.3, title, fontsize=20, horizontalalignment='center', verticalalignment='top',
    #                  color="white")
    #
    #     fig.savefig(output, bbox_inches='tight')

    @staticmethod
    def _plot_histogram(data, title="Fragment Hotspot Maps"):
            """
            initialise the matplotlib figure to output histograms

            :param data:
            :return:
            """
            colour_dict = {"acceptor": "r",
                           "donor": "b",
                           "apolar": "y",
                           "negative": "m",
                           "positive": "c"}
            plt.figure(1)
            for n, key in enumerate(data.keys()):
                j = int(len(data.keys()))
                plt.subplot(j, 1, (n + 1))
                hist, bin_edges = np.histogram(data[key], bins=range(0, 40), normed=True)
                Figures._histogram_settings(bin_edges)
                if n == 0:
                    plt.title(title)
                if n < j - 1:
                    plt.xticks([])
                if n == j - 1:
                    plt.xlabel("Fragment hotspot score")
                if n == int(round(j / 2) - 1):
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
        plt.xlim(min(bin_edges), max(bin_edges))
        plt.ylim(0, 0.35)
        plt.yticks([])

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

