"""
The :mod:`hotspots.hs_pharmacophore` module contains classes for the
conversion of Grid objects to pharmacophore models.

The main class of the :mod:`hotspots.hs_pharmacophore` module is:

- :class:`hotspots.hs_pharmacophore.PharmacophoreModel`


A Pharmacophore Model can be generated directly from a :class:`hotspots.result.Result` :

>>> from hotspots.calculation import Runner

>>> r = Runner()
>>> result = r.from_pdb("1hcl")
>>> result.get_pharmacophore_model(identifier="MyFirstPharmacophore")

The Pharmacophore Model can be used in Pharmit or CrossMiner

>>> result.pharmacophore.write("example.cm")   # CrossMiner
>>> result.pharmacophore.write("example.json")    # Pharmit

More information about CrossMiner is available:
    - Korb O, Kuhn B, hert J, Taylor N, Cole J, Groom C, Stahl M "Interactive and Versatile Navigation of Structural Databases" J Med Chem, 2016, 59(9):4257, [DOI: 10.1021/acs.jmedchem.5b01756]

More information about Pharmit is available:
    - Jocelyn Sunseri, David Ryan Koes; Pharmit: interactive exploration of chemical space, Nucleic Acids Research, Volume 44, Issue W1, 8 July 2016, Pages W442-W448 [DIO: 10.1093/nar/gkw287]

"""
from __future__ import print_function
import csv
import json
import re
import os
import tempfile
from os.path import basename, splitext, join, dirname

import numpy as np
from ccdc import io
from ccdc.descriptors import GeometricDescriptors
from ccdc.molecule import Atom, Molecule
from ccdc.pharmacophore import Pharmacophore
from ccdc.protein import Protein

from hotspots.grid_extension import Grid, Coordinates
from hotspots.hs_utilities import Helper
from hotspots.template_strings import pymol_arrow, pymol_imports, crossminer_features, pymol_labels
from hotspots.pdb_python_api import Query, PDB, PDBResult

from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, AllChem
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from tqdm import tqdm

from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import hdbscan



def tanimoto_dist(a, b):
    """
    calculate the tanimoto distance between two fingerprint arrays
    :param a:
    :param b:
    :return:
    """
    dotprod = np.dot(a, b)
    tc = dotprod / (np.sum(a) + np.sum(b) - dotprod)
    return 1.0 - tc


class _Ligand(object):
    """

    """
    def __init__(self, ccdc_mol, rdkit_mol, fingerprint, chem_id):
        self.ccdc_mol = ccdc_mol
        self.rdmol = rdkit_mol
        self.fingerprint = fingerprint
        self.chemical_id = chem_id

    @staticmethod
    def from_file(path):
        """


        :param path:
        :return:
        """
        ccdc_mol = io.MoleculeReader(path)[0]
        rdkit_mol = Chem.SDMolSupplier(path)[0]
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(rdkit_mol, 2)
        chem_id = ccdc_mol.identifier
        return _Ligand(ccdc_mol, rdkit_mol, fingerprint, chem_id)


class PharmacophoreModel(Helper):
    """
    A class to handle a Pharmacophore Model

    :param `hotspots.hs_pharmacophore.PharmacophoreModel.Settings` settings: Pharmacophore Model settings
    :param str identifier: Model identifier
    :param list features: list of  :class:hotspots.hs_pharmacophore._PharmacophoreFeatures
    :param `ccdc.protein.Protein` protein: a protein
    :param dict dic: key = grid identifier(interaction type), value = :class:`ccdc.utilities.Grid`

    """

    class Settings():
        """
        settings available for adjustment

        :param float feature_boundary_cutoff: The map score cutoff used to generate islands
        :param float max_hbond_dist: Furthest acceptable distance for a hydrogen bonding partner (from polar feature)
        :param float radius: Sphere radius
        :param bool vector_on: Include interaction vector
        :param float transparency: Set transparency of sphere
        :param bool excluded_volume:  If True, the CrossMiner pharmacophore will contain excluded volume spheres
        :param float binding_site_radius: Radius of search for binding site calculation, used for excluded volume
        """

        def __init__(self, feature_boundary_cutoff=5, max_hbond_dist=5, radius=1.0, vector_on=False, transparency=0.6,
                     excluded_volume=True, binding_site_radius=12):
            self.feature_boundary_cutoff = feature_boundary_cutoff
            self.max_hbond_dist = max_hbond_dist
            self.radius = radius  # set more intelligently
            self.transparency = transparency
            self.excluded_volume = excluded_volume
            self.binding_site_radius = binding_site_radius

            if vector_on:
                self.vector_on = 1
            else:
                self.vector_on = 0

    def __init__(self, settings, identifier=None, features=None, protein=None, dic=None):
        self.identifier = identifier
        self._features = features

        self.fname = None
        self.projected_dict = {"True": ["donor", "acceptor"], "False": ["negative", "positive", "apolar"]}
        self.protein = protein

        if settings == None:
            self.settings = self.Settings()
        else:
            self.settings = settings

        self.dic = dic

    def _usr_moment(self):
        """
        generates USR moments for shape analysis
        :return:
        """
        try:
            from hotspots.hs_utilities import _generate_usr_moment
            type_dic = {"apolar": [],
                        "acceptor": [],
                        "donor": [],
                        "positive": [],
                        "negative": []
                        }

            for feat in self.features:
                type_dic[feat.feature_type].append(feat.feature_coordinates)

            coords_list = [np.array(v) for v in type_dic.values() if len(v) != 0]
            return _generate_usr_moment(fcoords_list=coords_list)
        except ImportError:
            print("To use this feature you must have USR installed")

    @property
    def features(self):
        """
        Interaction features of the Pharmacophore Model

        :return: list of :class:`hotspots.hs_pharmacophore._PharmacophoreFeature instance
        """
        return self._features

    def _comparision_dict(self, feature_threshold=0):
        """
        converts pharmacophore into comparision dictionary

        :return: dic
        """
        d = {}
        self.rank_features(max_features=20, feature_threshold=feature_threshold)
        rank_dic = {"apolar": 0, "donor": 0, "acceptor": 0, "negative":0, "positive":0}
        for feat in self._features:
            d.update({feat.score_value: [feat.feature_type, feat.feature_coordinates, rank_dic[feat.feature_type]]})
            rank_dic[feat.feature_type] += 1

        return d

    def rank_features(self, max_features=4, feature_threshold=0, force_apolar=True):
        """
        orders features by score

        :param int max_features: maximum number of features returned
        :param float feature_threshold: only features above this value are considered
        :param force_apolar: ensures at least one point is apolar
        :return: list of features


        >>> from hotspots.hs_io import HotspotReader

        >>> result = HotspotReader("out.zip").read()
        >>> model = result.get_pharmacophore_model()
        >>> print(len(model.features))
        38
        >>> model.rank_features(max_features=5)
        >>> print(len(model.features))
        5

        """
        if force_apolar:
            apolar_dict = {feat.score_value: feat for feat in self.features if feat.feature_type == "apolar"}
            top_apolar = sorted(apolar_dict.items(), key=lambda x: x[0], reverse=True)[0]
            apolar = apolar_dict[top_apolar[0]]

            score_dic = {feat.score_value: feat for feat in self.features if feat.feature_type != "apolar"}
            sorted_scores = sorted(score_dic.items(), key=lambda x: x[0], reverse=True)
            max_features -= 1
            if max_features > len(self.features):
                ordered_features = [feat for score, feat in sorted_scores if score > feature_threshold]
            else:
                ordered_features = [feat for score, feat in sorted_scores if score > feature_threshold][:max_features]

            ordered_features.append(apolar)

        else:
            score_dic = {feat.score_value: feat for feat in self.features}
            sorted_scores = sorted(score_dic.items(), key=lambda x: x[0], reverse=True)
            if max_features > len(self.features):
                ordered_features = [feat for score, feat in sorted_scores if score > feature_threshold]
            else:
                ordered_features = [feat for score, feat in sorted_scores if score > feature_threshold][:max_features]

        self._features = ordered_features

    def _get_pymol_pharmacophore(self, lfile):
        """
        create the pymol visualisation for a pharmacophore
        """
        pymol_out = """
cluster_dict = {{"{0}":[], "{0}_arrows":[]}}
""".format(self.identifier)
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

        pymol_out += pymol_labels(fname=lfile,
                                  objname="label_threshold_{}".format(self.identifier))

        pymol_out += """\ncmd.group('Pharmacophore_{0}', members= 'label_threshold_{0}')\n""".format(self.identifier)

        return pymol_out

    def _as_grid(self, feature_type=None, tolerance=2):
        """
        returns _features as grid
        """
        if feature_type == None:
            filtered_features = self._features

        else:
            filtered_features = [feat for feat in self._features if feat.feature_type == feature_type]

        x = [feat.feature_coordinates.x for feat in filtered_features]
        y = [feat.feature_coordinates.y for feat in filtered_features]
        z = [feat.feature_coordinates.z for feat in filtered_features]

        origin = [min(x) - tolerance, min(y) - tolerance, min(z) - tolerance]
        far_corner = [max(x) + tolerance, max(y) + tolerance, max(z) + tolerance]
        grd = Grid(origin=origin, far_corner=far_corner, spacing=0.5, default=0, _grid=None)

        for feat in filtered_features:
            grd.set_sphere(point=feat.feature_coordinates, radius=self.settings.radius, value=1, scaling='None')
        return grd

    def _get_binding_site_residues(self):
        """
        return a protein binding site
        """
        centroid = [feat.feature_coordinates for feat in self._features if feat.feature_type == "apolar"][0]
        prot = self.protein

        bs = prot.BindingSiteFromPoint(protein=self.protein,
                                       origin=centroid,
                                       distance=self.settings.binding_site_radius)

        bs_residues = [str(r.identifier) for r in bs.residues]
        protein_residues = [str(p.identifier) for p in prot.residues]
        deletes = list(set(protein_residues) - set(bs_residues))
        for delete in deletes:
            prot.remove_residue(delete)

        return prot

    def _get_crossminer_pharmacophore(self):
        """
        convert a PharmacophoreModel into a crossminer pharmacophore
        """
        # TODO: UPDATE WITH CHARGED FEATURES
        supported_features = {"acceptor_projected": "acceptor",
                              "donor_projected": "donor",
                              "ring": "apolar"}
        try:
            Pharmacophore.read_feature_definitions()
        except:
            raise ImportError("Crossminer is only available to CSD-Discovery")

        feature_definitions = {supported_features[fd.identifier]: fd for fd in
                               Pharmacophore.feature_definitions.values()
                               if fd.identifier in supported_features.keys()}

        model_features = []
        for feat in self._features:
            if feat.feature_type == "negative" or feat.feature_type == "positive":
                print("Charged feature not currently supported in CrossMiner: Its on the TODO list")

            else:
                sphere = GeometricDescriptors.Sphere(feat.feature_coordinates, self.settings.radius)

                if feat.projected_coordinates:
                    projected = GeometricDescriptors.Sphere(feat.projected_coordinates, self.settings.radius)
                    p = Pharmacophore.Feature(feature_definitions[feat.feature_type], *[sphere, projected])

                else:
                    p = Pharmacophore.Feature(feature_definitions[feat.feature_type], sphere)

                model_features.append(p)

        if self.settings.excluded_volume:
            if not self.protein:
                print("Pharmacophore Model must have protein to calculate excluded volume")
            else:
                bs = self._get_binding_site_residues()

                for residue in bs.residues:
                    mol = None
                    mol = Molecule(identifier="temp_residue")

                    # for a in residue.backbone_atoms:
                    #     ev = Pharmacophore.ExcludedVolume(GeometricDescriptors.Sphere(a.coordinates, 2))
                    #     model_features.append(ev)
                    for a in residue.backbone_atoms:
                        mol.add_atom(a)

                    centre = mol.centre_of_geometry()
                    ev = Pharmacophore.ExcludedVolume(GeometricDescriptors.Sphere(centre, 2))
                    model_features.append(ev)

        return Pharmacophore.Query(model_features)

    @staticmethod
    def _run_query(accession_id):
        """
        search the PDB for entries which have the same UniProt code, contain ligands and are
        within a resolution cutoff

        :param str accession_id: the accession ID for the target
        :return list: list of :class:`pdb_python_api.PDBResult` instances
        """
        # create query
        q = Query()
        q.add_term(query_type="UpAccessionIdQuery",
                   query_parameters={'accessionIdList': '{}'.format(accession_id)})

        q.add_term(query_type="NoLigandQuery",
                   query_parameters={'haveLigands': 'yes'},
                   conjunction='and')

        q.add_term(query_type="ResolutionQuery",
                   query_parameters={'refine.ls_d_res_high.comparator': 'between',
                                     'refine.ls_d_res_high.min': '0',
                                     'refine.ls_d_res_high.max': '2.5'},
                   conjunction='and')

        return PDB.search(q.query)

    @staticmethod
    def _get_ligands(results):
        """
        smiles to RDKit molecule

        :param list: list of :class:`pdb_python_api.PDBResult` instances
        :return: RDKit molecules
        """
        ligs = []
        uniques = []
        for entry in results:
            for l in entry.filtered_ligands:
                try:
                    l.rdmol = Chem.MolFromSmiles(l.smiles)
                    AllChem.Compute2DCoords(l.rdmol)
                    l.rdmol.SetProp("_Name", str(entry.identifier + "/" + l.chemical_id))
                    # l.fingerprint = MACCSkeys.GenMACCSKeys(l.rdmol)
                    l.fingerprint = AllChem.GetMorganFingerprintAsBitVect(l.rdmol, 2)
                    if l.chemical_id not in uniques:
                        ligs.append(l)
                        uniques.append(l.chemical_id)
                except:
                    continue
        return ligs

    @staticmethod
    def _cluster_ligands(ligands, t):
        """

        :return:
        """
        def fingerprint_array(ligands):
            X =[]
            for l in ligands:
                arr = np.zeros((0,))
                DataStructs.ConvertToNumpyArray(l.fingerprint, arr)
                X.append(arr)
            return X

        cluster_dic = {}

        # generate fingerprint array
        X = fingerprint_array(ligands)
        if len(X) < 2:
            X = fingerprint_array(ligands)
            if len(X) < 2:
                raise ValueError("Fingerprint array must contain more than 1 entry")

        # dimensionality reduction
        tsne_X = TSNE(n_components=2, metric=tanimoto_dist).fit_transform(np.array(X, dtype=np.float32))

        # clustering
        cluster_tsne = hdbscan.HDBSCAN(min_cluster_size=2, gen_min_span_tree=True)
        cluster_tsne.fit(tsne_X)

        for i, label in enumerate(cluster_tsne.labels_):
            if label == -1:
                continue
            else:
                if label in cluster_dic:
                    cluster_dic[label].append(ligands[i])
                else:
                    cluster_dic.update({label: [ligands[i]]})

        x = [tsne_X.T[0][j] for j, l in enumerate(cluster_tsne.labels_) if l != -1]
        y = [tsne_X.T[1][j] for j, l in enumerate(cluster_tsne.labels_) if l != -1]
        hue = [l for j, l in enumerate(cluster_tsne.labels_) if l != -1]

        seen = [-1]
        sx = []
        sy = []
        for k, l in enumerate(cluster_tsne.labels_):
            if l in seen:
                continue
            else:
                sx.append(tsne_X.T[0][k])
                sy.append(tsne_X.T[1][k])
                seen.append(l)

        plt.scatter(x, y, c=hue, cmap='RdBu', alpha=0.7)
        plt.scatter(sx, sy, c="black", marker="x")

        plt.title("{} clusters".format(t))
        plt.savefig("{}.png".format(t))
        plt.close()

        if len(cluster_dic) == 0:
            print("NO CLUSTERS FOUND")
            try:
                unique = {}
                for l in ligands:
                    hetid = l.chemical_id.split("_")[0]
                    if not hetid in unique:
                        unique.update({hetid: l})

                ligands = unique.values()
                cluster_dic = {i: [ligands[i]] for i in range(0, len(ligands))}
            except:
                cluster_dic = {i: [ligands[i]] for i in range(0, len(ligands))}

        return cluster_dic

    @staticmethod
    def _align_proteins(reference, reference_chain, targets):
        """
        align proteins by chain

        :param `ccdc.protein.Protein` reference: align to me
        :param str reference_chain: align to this chain
        :param list targets: list of `ccdc.protein.Protein`
        :return tup: list(:class:`ccdc.protein.Protein`) and list (:classa`ccdc.molecule.Molecule`)
        """
        print("Aligning proteins to {}, chain {}...".format(reference.identifier, reference_chain))
        aligned_prots = []
        aligned_ligands = []

        reference = Protein.from_file(reference.fname)
        reference.add_hydrogens()

        for t in tqdm(targets):
            prot = Protein.from_file(t.fname)
            prot.detect_ligand_bonds()
            prot.add_hydrogens()
            chain = None
            for l in prot.ligands:
                if str(t.clustered_ligand) == str(l.identifier.split(":")[1][0:3]):
                    try:
                        bs = Protein.BindingSiteFromMolecule(protein=prot,
                                                             molecule=l,
                                                             distance=6)
                        chain = bs.residues[0].identifier.split(":")[0]
                    except:
                        chain = None

            if chain is None:
                print("\n        {} failed! No chain detected".format(t.identifier))
                continue
            try:
                binding_site_superposition = Protein.ChainSuperposition()
                (bs_rmsd, bs_transformation) = binding_site_superposition.superpose(reference[reference_chain],
                                                                                    prot[chain])
                aligned_prots.append(prot)
                for lig in prot.ligands:
                    if str(t.clustered_ligand) == str(lig.identifier.split(":")[1][0:3]):
                        if chain == str(lig.identifier.split(":")[0]):
                            aligned_ligands.append(lig)
            except:
                print("\n        {} failed!".format(t.identifier))
                continue

        return aligned_prots, aligned_ligands

    def write(self, fname):
        """
        writes out pharmacophore. Supported formats:

        - ".cm" (*CrossMiner*),
        - ".json" (`Pharmit <http://pharmit.csb.pitt.edu/search.html/>`_),
        - ".py" (*PyMOL*),
        - ".csv",
        - ".mol2"

        :param str fname: path to output file
        """
        extension = splitext(fname)[1]

        if extension == ".cm":

            print("WARNING! Charged features not currently supported in CrossMiner!")
            pharmacophore = self._get_crossminer_pharmacophore()
            pharmacophore.write(fname)

        elif extension == ".csv":
            with open(fname, "wb") as csv_file:
                csv_writer = csv.writer(csv_file, delimiter=",")
                line = 'Identifier, Feature_type, x, y, z, score, ' \
                       'projected_x, projected_y, projected_z, ' \
                       'vector_x, vector_y, vector_z'

                for feature in self._features:
                    line += "{0},{1},{2},{3},{4},{5}".format(self.identifier,
                                                             feature.feature_type,
                                                             feature.feature_coordinates.x,
                                                             feature.feature_coordinates.y,
                                                             feature.feature_coordinates.z,
                                                             feature.score_value
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
                lfile = "label_threshold_{}.mol2".format(self.identifier)

                pymol_out = pymol_imports()
                pymol_out += pymol_arrow()
                lines = self._get_pymol_pharmacophore(lfile)
                pymol_out += lines
                pymol_file.write(pymol_out)

            label = self.get_label(self)
            with io.MoleculeWriter(join(dirname(fname), lfile)) as writer:
                writer.write(label)

        elif extension == ".json":
            with open(fname, "w") as pharmit_file:
                pts = []
                interaction_dic = {'apolar': 'Hydrophobic',
                                   'donor': 'HydrogenDonor',
                                   'acceptor': 'HydrogenAcceptor',
                                   'negative': 'NegativeIon',
                                   'positive': 'PositiveIon'
                                   }

                for feat in self._features:
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
            mol = Molecule(identifier="pharmacophore_model")
            atom_dic = {"apolar": 'C',
                        "donor": 'N',
                        "acceptor": 'O',
                        "negative": 'S',
                        "positve": 'H'}

            pseudo_atms = [Atom(atomic_symbol=atom_dic[feat.feature_type],
                                atomic_number=14,
                                coordinates=feat.feature_coordinates,
                                label=str(feat.score_value))
                           for feat in self.features]

            for a in pseudo_atms:
                mol.add_atom(a)

            with io.MoleculeWriter(fname) as w:
                w.write(mol)

        elif extension == ".grd":
            g = self._as_grid()
            g.write(fname)

        else:
            raise TypeError("""""{}" output file type is not currently supported.""".format(extension))

    @staticmethod
    def from_hotspot(result, identifier="id_01", threshold=5, min_island_size=5, settings=None):
        """
        creates a pharmacophore model from a Fragment Hotspot Map result

        (included for completeness, equivalent to `hotspots.result.Result.get_pharmacophore()`)

        :param `hotspots.result.Result` result: a Fragment Hotspot Maps result (or equivalent)
        :param str identifier: Pharmacophore Model identifier
        :param float threshold: values above this value
        :param `hotspots.hs_pharmacophore.PharmacophoreModel.Settings` settings: settings

        :return: :class:`hotspots.hs_pharmacophore.PharmacophoreModel`

        >>> from hotspots.calculation import Runner
        >>> from hotspots.hs_pharmacophore import PharmacophoreModel

        >>> r = Runner()
        >>> result = r.from_pdb("1hcl")
        >>> model = PharmacophoreModel(result, identifier="pharmacophore")

        """

        if not settings:
            settings = PharmacophoreModel.Settings()

        feature_list = [_PharmacophoreFeature.from_hotspot(island, probe, result.protein, settings)
                        for probe, g in result.super_grids.items()
                        for island in g.islands(threshold) if island.count_grid() >= min_island_size]

        return PharmacophoreModel(settings,
                                  identifier=identifier,
                                  features=feature_list,
                                  protein=result.protein)

    @staticmethod
    def _from_file(fname, protein=None, identifier=None, settings=None):
        """
        creates a pharmacophore model from file (only .cm supported)

        :param str fname: path to CrossMiner file
        :param `ccdc.protein.Protein` protein: protein
        :param str identifier: identifier for the Pharmacophore Model
        :param `hotspots.hs_pharmacophore.PharmacophoreModel.Settings` settings: settings

        :return: :class:`hotspots.hs_pharmacophore.PharmacophoreModel`

        """
        # self, projected, feature_type, feature_coordinates, projected_coordinates, score_value, vector,
        #                  settings)
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
        if not settings:
            settings = PharmacophoreModel.Settings()

        if identifier is None:
            identifier = basename(fname).split(".")[0]

        # get feature information
        query = Pharmacophore.Query.from_file(fname)

        features = []
        for f in query.features:
            print(f.spheres[0].centre.point())
            centre = f.spheres[0].centre.point()
            if len(f.spheres) > 1:
                projected = True
                p = f.spheres[1].centre.point()
                features.append(_PharmacophoreFeature(projected=projected,
                                                      feature_type=cm_feature_dict[f.identifier],
                                                      feature_coordinates=Coordinates(float(centre[0]),
                                                                                      float(centre[1]),
                                                                                      float(centre[2])
                                                                                      ),
                                                      projected_coordinates=Coordinates(float(p[0]),
                                                                                        float(p[1]),
                                                                                        float(p[2])
                                                                                        ),
                                                      score_value=1,
                                                      vector=None,
                                                      settings=settings)
                                )
            else:
                projected = False
                features.append(_PharmacophoreFeature(projected=projected,
                                                      feature_type=cm_feature_dict[f.identifier],
                                                      feature_coordinates=Coordinates(float(centre[0]),
                                                                                      float(centre[1]),
                                                                                      float(centre[2])
                                                                                      ),
                                                      projected_coordinates=None,
                                                      score_value=1,
                                                      vector=None,
                                                      settings=settings)
                                )

        # with open(fname) as f:
        #     file = f.read().split("FEATURE_LIBRARY_END")[1]
        #     lines = [l for l in file.split("""\r\n\r\n""") if l != ""]
        #     feature_list = [f for f in [_PharmacophoreFeature.from_crossminer(feature) for feature in lines] if
        #                     f != None]

        return PharmacophoreModel(settings,
                                  identifier=identifier,
                                  features=features,
                                  protein=protein)

    @staticmethod
    def from_ligands(ligands, identifier, protein=None, settings=None):
        """
        creates a Pharmacophore Model from a collection of overlaid ligands

        :param `ccdc,molecule.Molecule` ligands: ligands from which the Model is created
        :param str identifier: identifier for the Pharmacophore Model
        :param `ccdc.protein.Protein` protein: target system that the model has been created for
        :param `hotspots.hs_pharmacophore.PharmacophoreModel.Settings` settings: Pharmacophore Model settings

        :return: :class:`hotspots.hs_pharmacophore.PharmacophoreModel`


        >>> from ccdc.io import MoleculeReader
        >>> from hotspots.hs_pharmacophore import PharmacophoreModel

        >>> mols = MoleculeReader("ligand_overlay_model.mol2")
        >>> model = PharmacophoreModel.from_ligands(mols, "ligand_overlay_pharmacophore")
        >>> # write to .json and search in pharmit
        >>> model.write("model.json")

        """
        cm_dic = crossminer_features()
        blank_grd = Grid.initalise_grid([a.coordinates for l in ligands for a in l.atoms])
        feature_dic = {"apolar": blank_grd.copy(),
                       "acceptor": blank_grd.copy(),
                       "donor": blank_grd.copy()}

        if not settings:
            settings = PharmacophoreModel.Settings()

        if isinstance(ligands[0], Molecule):
            temp = tempfile.mkdtemp()

            with io.MoleculeWriter(join(temp, "ligs.mol2")) as w:
                for l in ligands:
                    w.write(l)
            ligands = list(io.CrystalReader(join(temp, "ligs.mol2")))

        try:
            Pharmacophore.read_feature_definitions()
        except:
            raise ImportError("Crossminer is only available to CSD-Discovery")

        feature_definitions = [fd for fd in Pharmacophore.feature_definitions.values() if
                               fd.identifier != 'exit_vector' and
                               fd.identifier != 'heavy_atom' and
                               fd.identifier != 'hydrophobe' and
                               fd.identifier != 'fluorine' and
                               fd.identifier != 'bromine' and
                               fd.identifier != 'chlorine' and
                               fd.identifier != 'iodine' and
                               fd.identifier != 'halogen']

        for fd in feature_definitions:
            detected = [fd.detect_features(ligand) for ligand in ligands]
            all_feats = [f for l in detected for f in l]

            if not all_feats:
                continue

            for f in all_feats:
                feature_dic[cm_dic[fd.identifier]].set_sphere(f.spheres[0].centre, f.spheres[0].radius, 1)

        features = []
        for feat, feature_grd in feature_dic.items():
            peaks = feature_grd.get_peaks(min_distance=4, cutoff=1)
            for p in peaks:
                coords = Coordinates(p[0], p[1], p[2])
                projected_coordinates = None
                if feat == "donor" or feat == "acceptor":
                    if protein:
                        projected_coordinates = _PharmacophoreFeature.get_projected_coordinates(feat,
                                                                                                coords,
                                                                                                protein,
                                                                                                settings)
                features.append(_PharmacophoreFeature(projected=None,
                                                      feature_type=feat,
                                                      feature_coordinates=coords,
                                                      projected_coordinates=projected_coordinates,
                                                      score_value=feature_grd.value_at_coordinate(coords,
                                                                                                  position=False),
                                                      vector=None,
                                                      settings=settings
                                                      )
                                )

        return PharmacophoreModel(settings,
                                  identifier=identifier,
                                  features=features,
                                  protein=protein,
                                  dic=feature_dic)

    @staticmethod
    def _to_file(ligands, out_dir, fname="ligands.dat"):
        with open(os.path.join(out_dir, fname), "w") as f:
            for l in ligands:
                f.write("{},{},{}\n".format(l.structure_id, l.chemical_id, l.smiles))


    @staticmethod
    def _from_siena(pdb, ligand, mode, identifier, out_dir=None):
        """
        creates a Pharmacophore Model from a PDB code using a binding site search

        (Previously this was done through from_pdb however doing a binding site search is much cleaner)

        :param pdb_list:
        :return:
        """

        from hotspots.siena import Search
        searcher = Search()
        ensemble = searcher.create_ensemble(pdb_code=pdb,
                                            ligand=ligand,
                                            mode=mode)

        if out_dir is None:
            out_dir = tempfile.mkdtemp()

        ensemble.save(out_dir=out_dir)
        ligands = [_Ligand.from_file(path=os.path.join(out_dir, "ligands", f))
                   for f in os.listdir(os.path.join(out_dir, "ligands")) if f.split(".")[1] == "sdf"]

        cluster_dict = PharmacophoreModel._cluster_ligands(ligands=ligands, t=identifier)
        reps = [l[0].ccdc_mol for l in cluster_dict.values() if len(l) != 0]

        p = PharmacophoreModel.from_ligands(ligands=reps, identifier=identifier)
        p.representatives = reps
        return p

    @staticmethod
    def from_pdb(pdb_code, chain, representatives=None, identifier="LigandBasedPharmacophore"):
        """
        creates a Pharmacophore Model from a PDB code.

        This method is used for the creation of Ligand-Based pharmacophores. The PDB is searched for protein-ligand
        complexes of the same UniProt code as the input. These PDB's are align, the ligands are clustered and density
        of atom types a given point is assigned to a grid.

        :param str pdb_code: single PDB code from the target system
        :param str chain: chain of interest
        :param str out_dir: path to output directory
        :param representatives: path to .dat file containing previously clustered data (time saver)
        :param str identifier: identifier for the Pharmacophore Model

        :return: :class:`hotspots.hs_pharmacophore.PharmacophoreModel`


        >>> from hotspots.hs_pharmacophore import PharmacophoreModel
        >>> from hotspots.result import Results
        >>> from hotspots.hs_io import HotspotWriter
        >>> from ccdc.protein import Protein
        >>> from pdb_python_api import PDBResult


        >>> # get the PDB ligand-based Pharmacophore for CDK2
        >>> model = PharmacophoreModel.from_pdb("1hcl")

        >>> # the models grid data is stored as PharmacophoreModel.dic
        >>> # download the PDB file and create a Results
        >>> PDBResult("1hcl").download(<output_directory>)
        >>> result = Result(protein=Protein.from_file("<output_directory>/1hcl.pdb"), super_grids=model.dic)
        >>> with HotspotWriter("<output_directory>") as w:
        >>>     w.write(result)

        """
        temp = tempfile.mkdtemp()
        ref = PDBResult(pdb_code)
        ref.download(out_dir=temp, compressed=False)
        all_ligands = None
        if representatives:
            print("Reading representative PDB codes ...")
            reps = []
            f = open(representatives, "r")
            entries = f.read().splitlines()

            for entry in entries:
                pdb_code, hetid, smiles = entry.split(",")
                reps.append((pdb_code, hetid))

        else:
            accession_id = PDBResult(pdb_code).protein.sub_units[0].accession_id
            results = PharmacophoreModel._run_query(accession_id)
            # return all_ligands
            all_ligands = PharmacophoreModel._get_ligands(results)

            cluster_dict = PharmacophoreModel._cluster_ligands(ligands=all_ligands, t=identifier)
            reps = [l[0] for l in cluster_dict.values() if len(l) != 0]

        targets = []

        for rep in reps:
            try:
                r = PDBResult(identifier=rep.structure_id)
                r.clustered_ligand = rep.chemical_id

            except:
                try:
                    r = PDBResult(identifier=rep[0])
                    r.clustered_ligand = rep[1]
                except:
                    raise AttributeError

            r.download(out_dir=temp, compressed=False)
            targets.append(r)

        prots, ligands = PharmacophoreModel._align_proteins(reference=ref,
                                                            reference_chain=chain,
                                                            targets=targets)

        p = PharmacophoreModel.from_ligands(ligands=ligands, identifier=identifier)
        p.all_ligands = all_ligands
        p.representatives = reps
        p.aligned_ligands = ligands

        return p


class _PharmacophoreFeature(Helper):
    """
    A class to construct pharmacophoric models based upon fragment hotspot maps.

    :param projected:
    :param feature_type:
    :param feature_coordinates:
    :param projected_coordinates:
    :param score_value:
    :param vector:
    :param settings:

    """

    def __init__(self, projected, feature_type, feature_coordinates, projected_coordinates, score_value, vector,
                 settings):
        self._projected = projected
        self._feature_type = feature_type
        self._feature_coordinates = feature_coordinates
        self._projected_coordinates = projected_coordinates
        self._score_value = score_value
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
    def score_value(self):
        return self._score_value

    @property
    def vector(self):
        return self._vector

    @staticmethod
    def from_hotspot(grid, probe, protein, settings):
        """
        create a feature from a hotspot island

        :param `ccdc.utilities.Grid`grid: input grid
        :param str probe: probe type identifier
        :param `ccdc.protein.Protein` protein: target protein
        :param `hotspots.hs_pharmacophore.PharmacophoreModel.Settings` settings: settings
        :return: :class:`hotspots.hs_pharmacophore._PharmacophoreFeature`
        """

        def score_feature(grid, threshold=14, percentile=50):
            """
            returns
            :return:
            """
            return grid.grid_score(threshold=threshold, percentile=percentile)

        feature_type = probe
        if probe == "apolar":
            score, feature_coordinates = _PharmacophoreFeature.get_centroid(grid)
            #turn on for median score
            #score = score_feature(grid)
            projected = False
            projected_coordinates = None
            vector = None

        else:
            vector = None
            projected_coordinates = None
            score, feature_coordinates = _PharmacophoreFeature.get_maxima(grid)
            # turn on for median score
            #score = score_feature(grid)
            if probe == "donor" or probe == "acceptor":

                projected = True
                if protein:
                    projected_coordinates = _PharmacophoreFeature.get_projected_coordinates(feature_type,
                                                                                            feature_coordinates,
                                                                                            protein,
                                                                                            settings)
                else:
                    projected_coordinates = None
                    if projected_coordinates:
                        vector = _PharmacophoreFeature.get_vector(projected_coordinates, feature_coordinates)
            else:
                projected = False

        return _PharmacophoreFeature(projected, feature_type, feature_coordinates, projected_coordinates, score, vector,
                                     settings)

    @staticmethod
    def from_crossminer(feature_str):
        """
        create feature from CrossMiner string

        :param str feature_str: extracted string from CrossMiner file
        :return: :class:`hotspots.hs_pharmacophore._PharmacophoreFeature`
        """
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
        settings = PharmacophoreModel.Settings()

        feat = re.search(r"""PHARMACOPHORE_FEATURE (.+?)\r""", feature_str)
        if feat.group(1) == "excluded_volume":
            pass
        else:
            feature_type = cm_feature_dict[feat.group(1)]

            spher = re.findall("""PHARMACOPHORE_SPHERE (.+?)\r""", feature_str)
            print(spher)

            if len(spher) == 1:
                coords = spher[0].split(" ")
                feature_coordinates = Coordinates(float(coords[0]), float(coords[1]), float(coords[2]))
                projected = False
                projected_coordinates = None
                score = coords[3]

            elif len(spher) == 2:
                coords = spher[0].split(" ")
                proj = spher[1].split(" ")
                feature_coordinates = Coordinates(float(coords[0]), float(coords[1]), float(coords[2]))
                projected = True
                projected_coordinates = Coordinates(float(proj[0]), float(proj[1]), float(proj[2]))
                score = coords[3]

            else:
                raise IOError("feature format not recognised")

            return _PharmacophoreFeature(projected, feature_type, feature_coordinates, projected_coordinates, score,
                                         vector,
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
            dist = Helper.get_distance(atm.coordinates, feature_coordinates)
            if dist < settings.max_hbond_dist:
                if dist in near_atoms.keys():
                    near_atoms[dist].append(atm)
                else:
                    near_atoms.update({dist: [atm]})
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
        print(max_value)
        indices = grid.indices_at_value(max_value)

        if len(indices) == 1:
            coords = grid.indices_to_point(indices[0][0], indices[0][1], indices[0][2])
            return max_value, Coordinates(coords[0], coords[1], coords[2])
        else:
            coords = grid.indices_to_point(round(sum(i[0] for i in indices) / len(indices)),
                                           round(sum(j[1] for j in indices) / len(indices)),
                                           round(sum(k[2] for k in indices) / len(indices))
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
