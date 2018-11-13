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
The :mod:`fragment_hotspot_maps.pharmacophore` module contains classes for the
conversion of Grid objects to pharmacophore models.

The main classes of the :mod:`fragment_hotspot_maps.pharmacophore` module are:

- :class:`fragment_hotspot_maps.pharmacophore.PharmacophoreFeature`
- :class:`fragment_hotspot_maps.pharmacophore.PharmacophoreModel`

TO DO:
Could this be made into and extension of the ccdc.pharmacophore :mod:
"""
import collections
from os.path import basename, splitext
from template_strings import pymol_arrow, pymol_imports
import numpy as np
import csv
import json
import re

from ccdc import io
from ccdc.molecule import Atom, Molecule
from ccdc.descriptors import GeometricDescriptors
from ccdc.pharmacophore import Pharmacophore

from utilities import Utilities
from grid_extension import Grid

Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


class Settings():
    """
    settings for the PharmacophoreFeature class
    """

    def __init__(self):
        """
        feature_boundary_cutoff is the value of the island cutoff used. (boundaries addressed in extract hotspot)
        max_hbond_dist is the furtherest acceptable distance for a hydrogen bonding partner (from polar feature)
        """
        self.feature_boundary_cutoff = 5
        self.max_hbond_dist = 5
        self.radius = 1.0    # set more intelligently
        self.vector_on = 0
        self.transparency = 0.6
        self.excluded_volume = True
        self.binding_site_radius = 12


class PharmacophoreModel(object):
    """
    A class to wrap pharmacophore features and output in various formats
    """

    def __init__(self, settings, identifier=None, features=None, protein=None):
        """
        identifier is useful for displaying multiple models at once
        :param features:
        """
        self.identifier = identifier
        self._features = features

        self.fname = None
        self.projected_dict = {"True": ["donor", "acceptor"], "False": ["negative", "positive", "apolar"]}

        self.protein = protein
        if settings == None:
            self.settings = Settings()
        else:
            self.settings = settings

    @property
    def features(self):
        """list of `fragment_hotspot_maps.pharmacophore.PharmacophoreFeature` class instances"""
        return self._features

    def rank_polar_features(self, threshold=0, max_num_features=None):
        """
        rank polar features, sort feature list
        :return:
        """
        apolar = [feat for feat in self.features if feat.feature_type == "apolar"][0]
        score_dic = {feat.score: feat for feat in self.features if feat.feature_type != "apolar"}
        sorted_scores = sorted(score_dic.items(), key=lambda x: x[0], reverse=True)

        if max_num_features is None:
            ordered_features = [feat for score, feat in sorted_scores if score > threshold]
        else:
            ordered_features = [feat for score, feat in sorted_scores if score > threshold][:max_num_features]
        return ordered_features

    def get_pymol_pharmacophore(self):
        """

        :return:
        """
        pymol_out = """cluster_dict = {{"{0}":[], "{0}_arrows":[]}}""".format(self.identifier)
        sphere_dict = {'acceptor': '[COLOR, 1.00, 0.00, 0.00]',
                       'donor': '[COLOR, 0.00, 0.00, 1.00]',
                       'apolar': '[COLOR, 1.00, 1.000, 0.000]',
                       'surface': '[COLOR, 0.5, 0.5, 0.5]',
                       'positive': '[COLOR, 0.0, 1.0, 1.0]',
                       'negative': '[COLOR, 0.6, 0.1, 0.6]'
                       }
        colour_dict = {'acceptor': 'red blue',
                       'donor': 'blue red',
                       'apolar': 'yellow'}
        i = 0
        for feature in self.features:
            if feature.feature_type in self.projected_dict["True"] and feature.projected_coordinates is not None:
                i += 1
                arrow = 'cluster_dict["{7}_arrows"] += cgo_arrow([{0},{1},{2}], [{3},{4},{5}], color="{6}", name="Arrows_{7}_{8}")\n' \
                    .format(feature.feature_coordinates.x,
                            feature.feature_coordinates.y,
                            feature.feature_coordinates.z,
                            feature.projected_coordinates.x,
                            feature.projected_coordinates.y,
                            feature.projected_coordinates.z,
                            colour_dict[feature.feature_type],
                            self.identifier,
                            str(i))
            else:
                arrow = ''

            sphere = '{0} + [ALPHA, {1}] + [SPHERE, float({2}), float({3}), float({4}), float({5})]\n' \
                .format(sphere_dict[feature.feature_type],
                        feature.settings.transparency,
                        feature.feature_coordinates.x,
                        feature.feature_coordinates.y,
                        feature.feature_coordinates.z,
                        feature.settings.radius)

            pymol_out += '\ncluster_dict["{0}"] += {1}'.format(self.identifier, sphere)
            pymol_out += '\n{}'.format(arrow)

        pymol_out += '\ncmd.load_cgo(cluster_dict["{0}"], "Features_{0}", 1)' \
                     '\ncmd.load_cgo(cluster_dict["{0}_arrows"], "Arrows_{0}")'.format(self.identifier)
        pymol_out += '\ncmd.set("transparency", 0.2,"Features_{0}")' \
                     '\ncmd.group("Pharmacophore_{0}", members="Features_{0}")' \
                     '\ncmd.group("Pharmacophore_{0}", members="Arrows_{0}")\n'.format(self.identifier)
        return pymol_out

    def as_grid(self, feature_type=None, tolerance=2):
        """
        returns features as grid
        :param tolerance:
        :return:
        """
        if feature_type == None:
            filtered_features = self.features

        else:
            filtered_features = [feat for feat in self.features if feat.feature_type == feature_type]

        x = [feat.feature_coordinates.x for feat in filtered_features]
        y = [feat.feature_coordinates.y for feat in filtered_features]
        z = [feat.feature_coordinates.z for feat in filtered_features]

        origin = [min(x) - tolerance, min(y) - tolerance, min(z) - tolerance]
        far_corner = [max(x) + tolerance, max(y) + tolerance, max(z) + tolerance]
        grd = Grid(origin=origin, far_corner=far_corner, spacing=0.5, default=0, _grid=None)

        for feat in filtered_features:
            grd.set_sphere(point=feat.feature_coordinates,radius=self.settings.radius, value=1,scaling='None')
        return grd

    def _get_binding_site_residues(self):
        """

        :return:
        """
        centroid = [feat.feature_coordinates for feat in self.features if feat.feature_type == "apolar"][0]
        prot = self.protein

        bs = prot.BindingSiteFromPoint(protein=self.protein,
                                       origin=centroid,
                                       distance=self.settings.binding_site_radius)

        #print "bs", len(bs.residues)

        bs_residues = [str(r.identifier) for r in bs.residues]
        protein_residues = [str(p.identifier) for p in prot.residues]
        deletes = list(set(protein_residues) - set(bs_residues))
        for delete in deletes:
            prot.remove_residue(delete)

        return prot

    def get_crossminer_pharmacophore(self):
        """

        :return:
        """
        # TODO: UPDATE WITH CHARGED FEATURES
        supported_features = {"acceptor_projected": "acceptor",
                              "donor_projected": "donor",
                              "ring": "apolar"}
        try:
            Pharmacophore.read_feature_definitions()
        except:
            raise ImportError("CrossMiner must be installed")

        feature_definitions = {supported_features[fd.identifier]: fd for fd in Pharmacophore.feature_definitions.values()
                               if fd.identifier in supported_features.keys()}

        model_features = []
        for feat in self.features:
            if feat.feature_type == "negative" or feat.feature_type == "positive":
                print "Charged feature not currently supported in CrossMiner: Its on the TODO list"

            else:
                sphere = GeometricDescriptors.Sphere(feat.feature_coordinates, self.settings.radius)

                if feat.projected_coordinates:
                    projected = GeometricDescriptors.Sphere(feat.projected_coordinates, self.settings.radius)
                    p = Pharmacophore.Feature(feature_definitions[feat.feature_type], *[sphere, projected])

                else:
                    p = Pharmacophore.Feature(feature_definitions[feat.feature_type], sphere)

                model_features.append(p)

        if self.settings.excluded_volume:
            if self.protein == None:
                print "Pharmacophore Model must have protein to calculate excluded volume"
            else:
                bs = self._get_binding_site_residues()

                for residue in bs.residues:
                    mol = None
                    mol = Molecule(identifier = "temp_residue")

                    # for a in residue.backbone_atoms:
                    #     ev = Pharmacophore.ExcludedVolume(GeometricDescriptors.Sphere(a.coordinates, 2))
                    #     model_features.append(ev)
                    for a in residue.backbone_atoms:
                        mol.add_atom(a)

                    centre = mol.centre_of_geometry()
                    ev = Pharmacophore.ExcludedVolume(GeometricDescriptors.Sphere(centre, 2))
                    model_features.append(ev)
        #print len(model_features)
        return Pharmacophore.Query(model_features)

    def write(self, fname):
        """
        given a fname, will output Pharmacophore in detected format
        :param fname: str, extensions support: ".cm", ".py", ".json", ".csv", ".mol2"
        :return:
        """
        extension = splitext(fname)[1]

        if extension == ".cm":
            print "WARNING! Charged features not currently supported in CrossMiner!"
            pharmacophore = self.get_crossminer_pharmacophore()
            pharmacophore.write(fname)

        elif extension == ".csv":
            with open(fname, "wb") as csv_file:
                csv_writer = csv.writer(csv_file, delimiter=",")
                line = 'Identifier, Feature_type, x, y, z, score, ' \
                       'projected_x, projected_y, projected_z, ' \
                       'vector_x, vector_y, vector_z'

                for feature in self.features:
                    line += "{0},{1},{2},{3},{4},{5}".format(self.identifier,
                                                            feature.feature_type,
                                                            feature.feature_coordinates.x,
                                                            feature.feature_coordinates.y,
                                                            feature.feature_coordinates.z,
                                                            feature.score
                                                            )
                    if feature.projected_coordinates:
                        line += ",{0},{1},{2}".format(feature.projected_coordinates.x,
                                                      feature.projected_coordinates.y,
                                                      feature.projected_coordinates.z)
                    else:
                        line += ",0,0,0"

                    if feature.vector:
                        line += ",{0},{1},{2}".format(feature.vector.x,
                                                      feature.vector.y,
                                                      feature.vector.z)
                    else:
                        line += ",0,0,0"

                    l = line.split(",")
                    csv_writer.writerow(l)

        elif extension == ".py":
            with open(fname, "wb") as pymol_file:
                pymol_out = pymol_imports()
                pymol_out += pymol_arrow()
                lines = self.get_pymol_pharmacophore()
                pymol_out += lines
                pymol_file.write(pymol_out)

        elif extension == ".json":
            with open(fname, "w") as pharmit_file:
                pts = []
                interaction_dic = {'apolar': 'Hydrophobic',
                                   'donor': 'HydrogenDonor',
                                   'acceptor': 'HydrogenAcceptor',
                                   'negative': 'NegativeIon',
                                   'positive': 'PositiveIon'
                                   }

                for feat in self.features:
                    if feat.vector:
                        point = {"name": interaction_dic[feat.feature_type],
                                 "hasvec": True,
                                 "x": feat.feature_coordinates.x,
                                 "y": feat.feature_coordinates.y,
                                 "z": feat.feature_coordinates.z,
                                 "radius": feat.settings.radius,
                                 "enabled": True,
                                 "vector_on": feat.settings.vector_on,
                                 "svector": {"x": feat.vector.x,
                                             "y": feat.vector.y,
                                             "z": feat.vector.z},
                                 "minsize": "",
                                 "maxsize": "",
                                 "selected": False
                                 }
                    else:
                        point = {"name": interaction_dic[feat.feature_type],
                                 "hasvec": False,
                                 "x": feat.feature_coordinates.x,
                                 "y": feat.feature_coordinates.y,
                                 "z": feat.feature_coordinates.z,
                                 "radius": feat.settings.radius,
                                 "enabled": True,
                                 "vector_on": feat.settings.vector_on,
                                 "svector": {"x": 0,
                                             "y": 0,
                                             "z": 0},
                                 "minsize": "",
                                 "maxsize": "",
                                 "selected": False
                                 }
                    pts.append(point)
                pharmit_file.write(json.dumps({"points": pts}))

        elif extension == ".mol2":
            mol = Molecule(identifier = "pharmacophore_model")
            atom_dic = {"apolar": 'C',
                        "donor": 'N',
                        "acceptor": 'O',
                        "negative": 'S',
                        "positve": 'H'}

            pseudo_atms = [Atom(atomic_symbol=atom_dic[feat.feature_type],
                                atomic_number=14,
                                coordinates=feat.feature_coordinates,
                                label = str(feat.score))
                           for feat in self.features]

            for a in  pseudo_atms:
                mol.add_atom(a)

            with io.MoleculeWriter(fname) as w:
                w.write(mol)

        else:
            raise TypeError("""""{}" output file type is not currently supported.""".format(extension))

    @staticmethod
    def from_hotspot(protein, super_grids, identifier="id_01", cutoff=5, settings=None):
        """creates a pharmacophore model from hotspot results object"""
        settings = Settings()
        feature_list = [PharmacophoreFeature.from_hotspot(island, probe, protein, settings)
                        for probe, g in super_grids.items()
                        for island in g.islands(cutoff) if island.count_grid() >= 5]

        return PharmacophoreModel(settings, identifier=identifier, features=feature_list, protein=protein)

    @staticmethod
    def from_file(fname, protein=None, identifier=None, settings=None):
        """creates a pharmacophore model from file (only .cm supported) """
        if identifier == None:
            identifier = basename(fname).split(".")[0]

        with open(fname) as f:
            file = f.read().split("FEATURE_LIBRARY_END")[1]
            lines = [l for l in file.split("""\n\n""") if l != ""]
            feature_list = [f for f in [PharmacophoreFeature.from_crossminer(feature) for feature in lines] if f != None]


        return PharmacophoreModel(settings, identifier=identifier, features=feature_list, protein=protein)


class PharmacophoreFeature(Utilities):
    """

    A class to construct pharmacophoric models based upon fragment hotspot maps.
    This feature is designed to be used after fragment sized hotspots have been extracted.
    (Hotspot.extract_hotspots method)
    """

    def __init__(self, projected, feature_type, feature_coordinates, projected_coordinates, score, vector, settings):
        """
        new __init__ allows flexibility of pharmacophore generation
        :return:
        """

        self._projected = projected
        self._feature_type = feature_type
        self._feature_coordinates = feature_coordinates
        self._projected_coordinates = projected_coordinates
        self._score = score
        self._vector = vector

        self.settings = settings

        if feature_type == "donor" or "acceptor":
            self.settings.vector_on = 1

    @property
    def projected(self):
        return self._projected

    @property
    def feature_type(self):
        return self._feature_type

    @property
    def feature_coordinates(self):
        return self._feature_coordinates

    @property
    def projected_coordinates(self):
        return self._projected_coordinates

    @property
    def score(self):
        return self._score

    @property
    def vector(self):
        return self._vector


    @staticmethod
    def from_hotspot(grid, probe, protein, settings):
        """

        :param grid:
        :param probe:
        :param protein:
        """

        feature_type = probe
        if probe == "apolar":
            score, feature_coordinates = PharmacophoreFeature.get_centroid(grid)
            projected = False
            projected_coordinates = None
            vector = None

        else:
            vector = None
            projected_coordinates = None
            score, feature_coordinates = PharmacophoreFeature.get_maxima(grid)
            if probe == "donor" or probe == "acceptor":

                projected = True
                if protein:
                    projected_coordinates = PharmacophoreFeature.get_projected_coordinates(feature_type,
                                                                                           feature_coordinates,
                                                                                           protein,
                                                                                           settings)
                else:
                    projected_coordinates = None
                    if projected_coordinates:
                        vector = PharmacophoreFeature.get_vector(projected_coordinates, feature_coordinates)
            else:
                projected = False

        return PharmacophoreFeature(projected, feature_type, feature_coordinates, projected_coordinates, score, vector,
                                    settings)

    @staticmethod
    def from_crossminer(feature_str):
        """generates a pharmacophore model from a crossminer file"""
        cm_feature_dict = {"ring": "apolar",
                           "ring_planar_projected": "apolar",
                           "ring_non_planar": "apolar",
                           "hydrophobic": "apolar",
                           "donor_ch_projected": "donor",
                           "donor_projected": "donor",
                           "donor": "donor",
                           "acceptor_projected": "acceptor",
                           "acceptor": "acceptor",
                           "negative": "",
                           "positive": "",
                           }

        vector = None
        settings = Settings()

        feat = re.search(r"""PHARMACOPHORE_FEATURE (.+?)\n""", feature_str)
        if feat.group(1) == "excluded_volume":
            pass
        else:
            feature_type =  cm_feature_dict[feat.group(1)]

            spher = re.findall("""PHARMACOPHORE_SPHERE (.+?)\n""", feature_str)

            if len(spher) == 1:
                coords = spher[0].split(" ")
                feature_coordinates = Coordinates(float(coords[0]), float(coords[1]), float(coords[2]))
                projected = False
                projected_coordinates = None
                score = coords[3]

            elif len(spher) == 2:
                coords = spher[0].split(" ")
                proj =  spher[1].split(" ")
                feature_coordinates = Coordinates(float(coords[0]), float(coords[1]), float(coords[2]))
                projected = True
                projected_coordinates = Coordinates(float(proj[0]), float(proj[1]), float(proj[2]))
                score = coords[3]

            else:
                raise IOError("feature format not recognised")

            return PharmacophoreFeature(projected, feature_type, feature_coordinates, projected_coordinates, score, vector,
                                        settings)

    @staticmethod
    def get_vector(projected_coordinates, feature_coordinates):
        """
        generates vector to hydrogen bonding partner on a protein.
        :return: Coordinates (named tuple)
        """
        return Coordinates(projected_coordinates.x - feature_coordinates.x,
                           projected_coordinates.y - feature_coordinates.y,
                           projected_coordinates.z - feature_coordinates.z)

    @staticmethod
    def get_projected_coordinates(feature_type, feature_coordinates, protein, settings):
        """
        for a given polar feature, the nearest h-bonding partner on the protein is located.
        :param protein: a :class:`ccdc.protein.Protein` instance
        :return: feature_coordinates for hydrogen-bonding partner
        """
        if feature_type == 'donor':
            atms = [a for a in protein.atoms if a.is_acceptor]
        else:
            atms = [a for a in protein.atoms if a.is_donor]

        near_atoms = {}
        for atm in atms:
            dist = Utilities.get_distance(atm.coordinates, feature_coordinates)
            if dist < settings.max_hbond_dist:
                if dist in near_atoms.keys():
                    near_atoms[dist].append(atm)
                else:
                    near_atoms.update({dist:[atm]})
            else:
                continue
        if len(near_atoms.keys()) == 0:
            return None

        else:
            closest = sorted(near_atoms.keys())[0]
            select = near_atoms[closest][0]
            return select.coordinates

    @staticmethod
    def get_maxima(grid):
        """
        given a grid will return the max point
        :param grid:
        :return:
        """
        max_value = grid.extrema[1]
        indices = grid.indices_at_value(max_value)

        if len(indices) == 1:
            coords = grid.indices_to_point(indices[0][0], indices[0][1], indices[0][2])
            return max_value, Coordinates(coords[0], coords[1], coords[2])
        else:
            coords = grid.indices_to_point(round(sum(i[0] for i in indices)/len(indices)),
                                           round(sum(j[1] for j in indices)/len(indices)),
                                           round(sum(k[2] for k in indices)/len(indices))
                                           )
            return max_value, Coordinates(coords[0], coords[1], coords[2])

    @staticmethod
    def get_centroid(grid):
        """
        given a grid will return the centre of mass(mass ==HS score) and the score at that point
        :return: score, float, and coords, ccdc.molecule.Coordinate
        """

        weighted_x = 0
        weighted_y = 0
        weighted_z = 0
        total_mass = 0
        nx, ny, nz = grid.nsteps

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    grid_val = grid.value(i, j, k)
                    x, y, z = grid.indices_to_point(i, j, k)
                    weighted_x += grid_val * x
                    weighted_y += grid_val * y
                    weighted_z += grid_val * z
                    total_mass += grid_val

        coords = Coordinates(np.divide(weighted_x, total_mass),
                                      np.divide(weighted_y, total_mass),
                                      np.divide(weighted_z, total_mass)
                                      )
        score = grid.value_at_point(coords)

        return float(score), coords
