import os
import tempfile
import colorsys
from collections import OrderedDict

from pprint import pprint
from scipy.spatial import distance
import numpy as np

from ccdc.pharmacophore import Pharmacophore
from ccdc.io import csd_directory, CrystalReader, MoleculeWriter, MoleculeReader
from ccdc.molecule import Molecule, Atom, Coordinates

from ccdc.descriptors import GeometricDescriptors
from ccdc import utilities
from ccdc.utilities import PushDir

from hotspots.grid_extension import Grid
from hotspots.wrapper_pdb import PDBResult
from hotspots.wrapper_arpeggio import Arpeggio
from hotspots.protein_extension import Protein
from hotspots.wrapper_pymol import PyMOLCommands, PyMOLFile

with utilities._private_importer():
    import MotifPharmacophoreLib
    import SubstructureSearchLib
    import ChemistryLib
    import MathsLib
    import UtilitiesLib


def _coordinate_str(coord, digits=4):
    return f"{round(coord[0], digits)}_{round(coord[1], digits)}_{round(coord[2], digits)}"


def to_array(t):
    return np.array([float(t[0]), float(t[1]), float(t[2])])


def rgb_to_decimal(colour):
    return utilities.Colour(r=round(colour.r / 255, 2),
                            g=round(colour.g / 255, 2),
                            b=round(colour.b / 255, 2),
                            a=round(colour.a / 255, 2))


def decimal_to_rgb(colour):
    return utilities.Colour(r=round(colour.r * 255, 2),
                            g=round(colour.g * 255, 2),
                            b=round(colour.b * 255, 2),
                            a=round(colour.a * 255, 2))


def adjust_lightness(colour, percentage):
    r, g, b, a = rgb_to_decimal(colour)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l += (percentage / 100)

    if l > 1:
        print("MAX lightness")
        l = 1.0

    elif l < 0:
        print("MIN lightness")
        l = 0.0

    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return utilities.Colour(r=round(r, 2),
                            g=round(g, 2),
                            b=round(b, 2),
                            a=round(a, 2))


def _create_grids(pharmacophores):
    g = Grid.initalise_grid(coords=[s.centre for p in pharmacophores for f in p.features for s in f.spheres],
                            spacing=0.5)
    fds = {f.identifier for p in pharmacophores for f in p.features}
    return {fd: g.copy() for fd in fds}


def _features_to_grid(all_features, grid, radius=2, scaling='linear', deduplicate=True):
    # PREVENT DOUBLE COUNTING
    # for each feature definition, only add the FeaturePoint (point of corresponding atom) once
    # for example, C=O Oxygen will generate 2 acceptor_projected features, one for each lone pair
    # we only want to count the oxygen atom position once. (It can be unpacked to 2 features later)

    if deduplicate:
        spheres = {(feature.spheres[0].centre[0],
                    feature.spheres[0].centre[1],
                    feature.spheres[0].centre[2]) for feature in all_features}
    else:
        spheres = [(feature.spheres[0].centre[0],
                    feature.spheres[0].centre[1],
                    feature.spheres[0].centre[2]) for feature in all_features]

    print(f"deduplicate: {len(all_features)} to {len(spheres)}")
    for sphere in spheres:
        grid.set_sphere(point=sphere, value=1, radius=radius, scaling=scaling)

    return grid


def select_projections(features, peak, tolerance=1):
    # projected spheres
    projs = np.array([[f.point[0].centre[0],
                       f.point[0].centre[1],
                       f.point[0].centre[2]] for f in features])

    # distance matrix
    l = distance.cdist(peak, projs)
    # create mask
    x, y = np.nonzero(l <= tolerance)

    # return features above threshold
    return [features[i] for i in set(y)]


def closest_peak_index(peaks_array, feature, max_distance):
    f = np.array([[feature.spheres[0].centre[0],
                   feature.spheres[0].centre[1],
                   feature.spheres[0].centre[2]]])

    d = distance.cdist(f, peaks_array)
    index = np.argmin(d)

    if d[0][index] < max_distance:
        return index
    else:
        return None


def create_consensus(pharmacophores, cutoff=2, max_distance=2.0):
    """

    """
    new_features = []

    # initialise grids from all pharmacophores
    feature_point_grids = _create_grids(pharmacophores)

    # add point spheres to corresponding grids
    features_by_type = {k: [f for p in pharmacophores for f in p.features if f.identifier == k]
                        for k in feature_point_grids.keys()}

    for identifier, all_features in features_by_type.items():
        # add spheres to grid
        feature_point_grids[identifier] = _features_to_grid(all_features, feature_point_grids[identifier])

        # find peaks
        all_peaks = feature_point_grids[identifier].get_peaks(min_distance=2, cutoff=0)

        peak_objs = []
        for peak in all_peaks:
            value = feature_point_grids[identifier].value_at_point(peak)
            if value >= cutoff:
                peak_objs.append(GridPeak(point=peak,
                                          value=value,
                                          feature_def=Pharmacophore.feature_definitions[identifier]))

        peaks_array = np.array([p.point for p in peak_objs])
        # projections
        if len(peaks_array) > 0:
            if len(all_features[0].spheres) > 1:
                # assign features to closest peak
                for feature in all_features:
                    index = closest_peak_index(peaks_array, feature, max_distance)
                    if index is not None:
                        peak_objs[int(index)].features.append(feature)

                # create projections
                for j, peak in enumerate(peak_objs):
                    peak.create_projection_grid()
                    peak.find_projection_peaks()
                    feats = peak.create_new_features()
                    new_features.extend(feats)
            else:
                for j, peak in enumerate(peak_objs):
                    point = GeometricDescriptors.Sphere(centre=peak.point, radius=1)
                    feat = Pharmacophore.Feature(peak.feature_def, point)
                    feat.point = point
                    feat.score = peak.value
                    new_features.append(feat)

    return new_features, feature_point_grids


class GridPeak:
    class ProjectionPeak:
        def __init__(self, point, value):
            self.point = point
            self.value = value

    def __init__(self, point, value, feature_def):
        self.point = point
        self.value = value
        self.feature_def = feature_def

        self.features = []

        # make from projections
        self.projection_grid = None
        self.projection_peaks = []
        self.new_features = []

    def __len__(self):
        return len(self.features)

    def create_projection_grid(self):
        """
        create the projected spheres for the summary pharmacophore features

        :return: `hotspots.grid_extension.Grid`
        """
        # Do the features have projections ?
        if len([f for f in self.features if len(f.spheres) > 1]) == 0:
            return None
        else:
            g = Grid.initalise_grid(coords=[feature.spheres[1].centre for feature in self.features],
                                    spacing=0.25,
                                    padding=2)

            for f in self.features:
                if len(f.spheres) > 1:
                    px, py, pz = [f.spheres[1].centre[0],
                                  f.spheres[1].centre[1],
                                  f.spheres[1].centre[2]]
                    g.set_sphere(point=[px, py, pz], value=1, radius=2, scaling='linear')
            self.projection_grid = g
            return g

    def find_projection_peaks(self):
        try:
            # the projection score has to be > 50 % of the point scire
            self.projection_peaks = [
                GridPeak.ProjectionPeak(point=peak, value=self.projection_grid.value_at_point(peak))
                for peak in self.projection_grid.get_peaks(min_distance=2, cutoff=1)
                if self.projection_grid.value_at_point(peak) > (0.5 * self.value)]
        except AttributeError:
            self.projection_peaks = []

    def create_new_features(self):
        """
        create new features from points and projections over the cutoff value
        :return: list of (`hotspots.pharmacophore_extension.Feature`)
        """
        new_feats = []
        point = GeometricDescriptors.Sphere(centre=self.point, radius=1)
        if len(self.projection_peaks) == 0:
            feat = Pharmacophore.Feature(self.feature_def, point)
            feat.point = point
            feat.score = self.value
            new_feats.append(feat)

        else:
            for projection in self.projection_peaks:
                proj = GeometricDescriptors.Sphere(centre=projection.point, radius=1)
                feat = Pharmacophore.Feature(self.feature_def, point, proj)
                feat.point = point
                feat.score = self.value
                feat.projected = proj
                feat.projected_value = projection.value
                new_feats.append(feat)
        return new_feats

    def to_pymol_str(self):
        pymol_out = ""
        pymol_out += PyMOLCommands.sphere("peak_obj", (1, 1, 1, 1), coords=self.point, radius=1)
        pymol_out += PyMOLCommands.load_cgo("peak_obj", "peak")

        ligand_peak_features = []
        for i, feat in enumerate(self.features):
            pymol_out += feat.to_pymol_str(i=i)
            ligand_peak_features.append(f"{feat.identifier}_{i}")

        pymol_out += PyMOLCommands.group("ligand_peak_features", ligand_peak_features)

        return pymol_out


class PharmacophoreModel(Pharmacophore.Query):
    """
    Base Class for the representing a CrossMiner pharmacophore query. Used to
    generate more specific PharmacophoreModels
    """

    class Settings:
        def __init__(self):
            self.x = 1

    def __init__(self):
        super().__init__()
        self.cm_dir = os.path.dirname(os.path.dirname(csd_directory()))
        Pharmacophore.read_feature_definitions(directory=os.path.join(self.cm_dir,
                                                                      "CSD_CrossMiner/feature_definitions"))
        self.__feature_options = {k: v for k, v in Pharmacophore.feature_definitions.items()}
        assert len(self.__feature_options) > 1

        self.__feature_definitions = self.__feature_options

        self.tmp = tempfile.mkdtemp()
        self.__identifier = None
        self.__ligands = None
        self.__protein = None
        self.__detected_features = None
        self.__feature_point_grids = None

    def __len__(self):
        return len(self.__features)

    @property
    def identifier(self):
        return self.__identifier

    @identifier.setter
    def identifier(self, ident):
        self.__identifier = ident

    @property
    def ligands(self):
        return self.__ligands

    @ligands.setter
    def ligands(self, ligands):
        self.__ligands = ligands

    @property
    def protein(self):
        return self.__protein

    @protein.setter
    def protein(self, protein):
        self.__protein = protein

    @property
    def feature_point_grids(self):
        return self.__feature_point_grids

    @feature_point_grids.setter
    def feature_point_grids(self, dict):
        self.__feature_point_grids = dict

    @property
    def detected_features(self):
        return self.__detected_features

    @detected_features.setter
    def detected_features(self, features):
        self.__detected_features = features

    @property
    def feature_definitions(self):
        return self.__feature_definitions

    @feature_definitions.setter
    def feature_definitions(self, feature_types):
        # reset before choosing
        self.__feature_definitions = {k: v for k, v in self.__feature_options.items()
                                      if any([k == ft for ft in feature_types])}

    def top_features(self, num, point=True, projection=True):
        feature_by_score = {}
        for feature in self.detected_features:
            score = []
            if point:
                score.append(feature.score)
            elif projection:
                score.append(feature.projected_value)

            feature_by_score.update({feature: sum(score)})

        print(feature_by_score)

        return list(OrderedDict(sorted(feature_by_score.items(), key=lambda item: item[1], reverse=True)).keys())[:num]

    def features_to_pymol_strings(self, features):
        """
        creates the code to visualise `ccdc.pharmacophore.Pharmacophore.Features` in PyMOL

        :param list of (`ccdc.pharmacophore.Pharmacophore.Features`) features: features to be visualised
        :return: str python code
        """
        group_dic = {ident: [] for ident in self.feature_definitions.keys()}
        pymol_out = ''

        for i, feat in enumerate(features):
            pymol_out += feat.to_pymol_str(i=i)
            group_dic[feat.identifier].append(f"{feat.identifier}_{i}")

            if feat.projected_identifier and self.protein:
                resnum = feat.projected_identifier.split("/")[1]
                pymol_out += f'\ncmd.select("sele", "resi {resnum}")\ncmd.show("sticks", "sele")'

        for fd in self.feature_definitions.keys():
            pymol_out += PyMOLCommands.group(f"{fd}_pts", group_dic[fd])
        pymol_out += PyMOLCommands.group("ligand_pharmacophore", [f"{a}_pts" for a in group_dic.keys()])

        return pymol_out

    def pymol_visulisation(self, outdir=None, fname="pymol_file.py"):
        if not outdir:
            outdir = os.getcwd()

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        if self.ligands:
            with MoleculeWriter(os.path.join(outdir, "ligands.mol2")) as w:
                for ligand in self.ligands:
                    try:
                        w.write(ligand.molecule)
                    except AttributeError:
                        w.write(ligand)

        if self.protein:
            with MoleculeWriter(os.path.join(outdir, "protein.mol2")) as w:
                w.write(self.protein.molecule)

        self.pymol_out = PyMOLFile()

        if self.ligands:
            self.pymol_out.commands += PyMOLCommands.load("ligands.mol2", "ligands")

        if self.protein:
            self.pymol_out.commands += PyMOLCommands.load("protein.mol2", "protein")

        # write out point spheres and projection sphere and lines if applicable
        self.pymol_out.commands += self.features_to_pymol_strings(self.detected_features)

        if self.feature_point_grids:
            for identifier, g in self.feature_point_grids.items():
                g.write(os.path.join(outdir, f"{identifier}.grd"))

                point_colour = rgb_to_decimal(self.feature_definitions[identifier].colour)
                self.pymol_out.commands += PyMOLCommands.set_color(f"{identifier}_color", point_colour)
                self.pymol_out.commands += PyMOLCommands.load(f"{identifier}.grd", f"{identifier}_grid")
                self.pymol_out.commands += PyMOLCommands.isosurface(f"{identifier}_grid",
                                                                    f"surface_{identifier}",
                                                                    level=1,
                                                                    color=f"{identifier}_color")

            self.pymol_out.commands += PyMOLCommands.group("feature_grids", self.feature_point_grids.keys())
            self.pymol_out.commands += PyMOLCommands.group("feature_grids",
                                                           [f"surface_{a}" for a in self.feature_point_grids])

            min_value = 0
            surface_dic = {
                self.identifier: {'feature_grids': [f"surface_{g}" for g in self.feature_point_grids.keys()]}}

            surface_value_dic = {self.identifier: {"feature_grids": max([round(g.extrema[1], 1)
                                                                         for g in self.feature_point_grids.values()])}}

            self.pymol_out.commands += PyMOLCommands.isoslider(surface_dic, surface_value_dic)

        self.pymol_out.write(os.path.join(outdir, fname))

    # STATIC METHODS
    @staticmethod
    def _get_crystal(obj):
        """
        Convert a obj's writable by MoleculeWriter to a crystal

        :param `ccdc.molecule.Molecule` obj: molecule or protein
        :return: `ccdc.crystal.Crystal`
        """
        tmp = tempfile.mkdtemp()
        f = os.path.join(tmp, "obj.mol2")
        with MoleculeWriter(f) as w:
            w.write(obj)
        return CrystalReader(f)[0]


class ProteinPharmacophoreModel(PharmacophoreModel):
    """
    Very simple derived class. Just takes advantage of some of the visualisation code
    """

    def __init__(self):
        super().__init__()

    def detect_from_prot(self, prot):
        """
        Add features to pharmacophore model from a ligand

        :param `ccdc.crystal.Crystal` or `ccdc.molecule.Molecule` ligand: a ligand
        :return:
        """
        if isinstance(prot, Protein):
            prot = self._get_crystal(prot)

        self.protein = prot
        self.detected_features = []
        for fd in self.feature_definitions.values():
            detected_feats = fd.detect_features(prot)
            if len(detected_feats) != 0:
                for f in detected_feats:
                    self.detected_features.append(f)


class LigandPharmacophoreModel(PharmacophoreModel):
    """
    Derived Class for representing pharmacophores generated from ligand atom types. Atom
    typing is carried out using built-in CrossMiner SMARTS definitions.

        - Pharmacophore generated from PDB codes use the RCSB PDB RESTful web service
        (https://www.rcsb.org/pdb/software/rest.do)

        - Ligand ensembles points are generated by collecting points on a grid and peak picking

    """

    def __init__(self):
        super().__init__()

    # from PDB
    def _align_proteins(self):
        return 0

    def _cluster_ligands(self):
        return 0

    # def add_self_features(self):
    #     # detection and addition separated to enable the selection to be adjusted after the
    #     # inital detection
    #     for feature in self.features:
    #         print(f"adding {feature.identifier} feature")
    #         self.add_feature(feature)

    def detect_from_ligand(self, ligand):
        """
        Add features to pharmacophore model from a ligand

        :param `ccdc.crystal.Crystal` or `ccdc.molecule.Molecule` ligand: a ligand
        :return:
        """
        if isinstance(ligand, Molecule):
            ligand = self._get_crystal(ligand)

        self.ligands = [ligand]

        self.detected_features = []
        for fd in self.feature_definitions.values():
            detected_feats = fd.detect_features(ligand)
            if len(detected_feats) != 0:
                for f in detected_feats:
                    self.detected_features.append(f)

    def detect_from_pdb(self, pdb, hetid, chainid=None):
        """
        Add features to pharmacophore model from a ligand in the PDB

        :param str pdb: pdb code
        :param str hetid: chemical component code
        :param str chainid: chain identifier
        :return: None
        """
        PDBResult(pdb).download(out_dir=self.tmp)
        prot = Protein.from_file(os.path.join(self.tmp, f"{pdb}.pdb"))
        prot.add_hydrogens()

        if chainid:
            lig = [l for l in prot.ligands
                   if l.identifier.split(":")[1][:3] == hetid and
                   l.identifier.split(":")[0] == chainid][0]

        else:
            lig = [l for l in prot.ligands if l.identifier.split(":")[1][:3] == hetid][0]

        self.detect_from_ligand(ligand=lig)


class HotspotPharmacophoreModel(PharmacophoreModel):
    """
    Derived Class for representing pharmacophores generated from Fragment Hotspot Maps.
    """

    def __init__(self):
        super().__init__()

    def from_hotspot(self, hr, projections=True):
        interaction_dict = {"donor": ["acceptor_projected"],
                            "acceptor": ["donor_projected",
                                         "donor_ch_projected"]
                            # "apolar": ["ring_planar_projected"]
                            }

        hotspot_to_cm = {"projected": {"apolar": "ring_planar_projected",
                                       "donor": "donor_projected",
                                       "acceptor": "acceptor_projected"},
                         "non-projected": {"apolar": "ring",
                                           "donor": "None",
                                           "acceptor": "acceptor"},
                         }

        # get peaks
        features = []
        for p, g in hr.super_grids.items():
            # peak as a sphere
            all_peaks = g.get_peaks(min_distance=1, cutoff=5)
            for peak in all_peaks:
                point = GeometricDescriptors.Sphere(centre=peak, radius=1)
                score = g.value_at_point(peak)

                if p != "apolar" and projections:
                    # binding site from point (within 4/5 angstrom of peak)
                    print("get binding site from point")
                    binding_site = hr.protein.copy()
                    bs = Protein.BindingSiteFromPoint(hr.protein, peak, distance=6)
                    for r in ({r.identifier for r in binding_site.residues} - {r.identifier for r in bs.residues}):
                        binding_site.remove_residue(r)

                    # detect projected features
                    print("detect cavity features")
                    pm = ProteinPharmacophoreModel()
                    pm.feature_definitions = interaction_dict[p]
                    pm.detect_from_prot(binding_site)
                    feats = pm.detected_features

                    # This returns multiple: ATM the user will then select which one (semi-automated)
                    # TODO: implement method to pick the best projection
                    projs = select_projections(feats, np.array([peak]), tolerance=4)

                    n_projs = len(projs)
                else:
                    n_projs = 0

                if n_projs == 0:
                    # no projections
                    if p == "donor":
                        print("Need to implement new CM feature def here, skipping for now")
                    else:
                        f = Pharmacophore.Feature(self.feature_definitions[hotspot_to_cm["non-projected"][p]],
                                                  point)
                        f.point = point
                        f.score = score
                        features.append(f)

                else:
                    for proj in projs:
                        centre = (proj.spheres[0].centre[0],
                                  proj.spheres[0].centre[1],
                                  proj.spheres[0].centre[2])

                        s = GeometricDescriptors.Sphere(centre=centre, radius=1)
                        f = Pharmacophore.Feature(self.feature_definitions[hotspot_to_cm["projected"][p]],
                                                  point,
                                                  s)
                        f.point = point
                        f.projected = s
                        f.score = score
                        features.append(f)

        self.detected_features = features


class InteractionPharmacophoreModel(PharmacophoreModel):
    """
    Derived Class for representing pharmacophores generated from Arpeggio.

    This approach uses detected protein-ligand interactions to generate pharmacophore features
    """

    def __init__(self):
        super().__init__()

    def detect_interactions(self, bs, ligand):

        interaction_partners = {"donor_projected": ["acceptor_projected", "ring"],
                                "donor_ch_projected": ["acceptor_projected", "ring"],
                                "acceptor_projected": ["donor_projected", "donor_ch_projected"],
                                "ring": ["ring", "donor_projected", "donor_ch_projected", "hydrophobe"],
                                "hydrophobe": ["ring"]}

        self.feature_definitions = ["acceptor_projected",
                                    "donor_projected",
                                    "donor_ch_projected",
                                    "ring",
                                    "hydrophobe"]

        if isinstance(ligand, Molecule):
            ligand = self._get_crystal(ligand)

        for fd in self.feature_definitions.values():
            detected_feats = fd.detect_features(ligand)
            print(fd.identifier, len(detected_feats))

        return 1

    def detect_from_arpeggio(self, protein_path, hetid, chain):
        """
        creates a pharmacophore from protein-ligand interactions
        TODO: This could be cleaner but for time reasons this is good enough. For example the SMARTS
              definitions between Arpeggio and Crossminer are not identical. Also, the SMARTS "grouping"
              are subtly different.

        1. Create atomtypes using Crossminer
        2. Detect bonds using Arepeggio
        3. For features in ligand, if bond, create a feature

        :param protein:
        :param hetid:
        :param chainid:
        :return:
        """
        # Arpeggio needs the protein in file, read the protein to get some information
        protein = Protein.from_file(protein_path)
        pdb_code = protein.identifier
        # assertion 1: protein must be protonated for CrossMiner
        # assert("H" in {atom.atomic_symbol for atom in protein.atoms[:50]})

        # assertion 2: protein must contain the ligand of interest
        assert (len([l for l in protein.ligands if l.identifier.split(":")[0] == chain and
                     l.identifier.split(":")[1][:3] == hetid]) == 1)

        lig = [l for l in protein.ligands if l.identifier.split(":")[0] == chain and
               l.identifier.split(":")[1][:3] == hetid][0]

        # CrossMiner needs a `ccdc.crystal.Crystal`
        crystal_ligand = self._get_crystal(lig)

        # Run Arpeggio
        arpeggio = Arpeggio(pdb_code, hetid, chain, protein_path)
        arpeggio.run()
        atom_features, ring_features = arpeggio.create_feature_list()
        interaction_features = []
        interaction_features.extend(atom_features)
        interaction_features.extend(ring_features)

        # CrossMiner atom-typing
        new_features = []
        ipoints = np.array([to_array(interaction.point) for interaction in interaction_features])
        for identifier, fd in self.feature_definitions.items():
            feats = fd.detect_features(crystal_ligand)

            fpoints = np.array([to_array(feat.point[0].centre) for feat in feats])

            # find the overlapping points
            ipoint_index, fpoint_index = np.where(distance.cdist(ipoints, fpoints) < 0.01)

            for i, f in zip(ipoint_index, fpoint_index):
                if identifier is "ring":
                    fd = Pharmacophore.feature_definitions["ring_planar_projected"]

                # sphere from CrossMiner
                point = GeometricDescriptors.Sphere(fpoints[f], 1)
                # sphere from Arpeggio
                projected = GeometricDescriptors.Sphere(to_array(interaction_features[i].projected), 1)
                print(point, projected)
                new = Pharmacophore.Feature(fd, point, projected)
                new.point = point
                new.point_identifier = interaction_features[i].point_identifier
                new.projected = projected
                new.projected_identifier = interaction_features[i].projected_identifier
                print(new.projected_identifier)

                new_features.append(new)

        print(len(new_features))

        self.detected_features = new_features
        self.protein = self._get_crystal(protein)
        self.ligands = [crystal_ligand]


class FeatureDefinition(Pharmacophore.FeatureDefinition):
    def __init__(self, _feature_def=None):
        super().__init__(_feature_def=_feature_def)

    def detect_features(self, crystal):
        _csv = ChemistryLib.CrystalStructureView_instantiate(crystal._crystal)
        _ssc = MotifPharmacophoreLib.MotifPharmacophoreSearchStructureCreator()
        _ssc.register_components_from_feature_definitions((self._feature_def,))
        _mss = _ssc.search_structure(_csv)
        _ded = self._feature_def.feature_deducer()
        _feats = _ded.generate_motif_features(_mss, self._feature_def.component_label())
        features = []
        for i in range(_feats.end_index()):
            code, feat_pts = _feats.at(i)
            for fp in feat_pts:
                pts = fp.points()
                spheres = tuple(GeometricDescriptors.Sphere((p[0], p[1], p[2]), 1.0) for p in pts)
                # Skip duplicates
                dup = False
                if len(spheres) == 2:
                    for f in features:
                        if len(f.spheres) == 2:
                            if (
                                    f.spheres[0] == spheres[1] and
                                    f.spheres[1] == spheres[0]
                            ):
                                dup = True
                                break
                if not dup:
                    feat = Pharmacophore.Feature(self._clone(), *spheres, crystal=crystal)
                    feat.label = f"{crystal.identifier}/{self.identifier}/{i}"
                    features.append(feat)
        return tuple(features)


Pharmacophore.FeatureDefinition = FeatureDefinition


class Feature(Pharmacophore.Feature):
    def __init__(self, feature_definition, *spheres, crystal=None):
        super().__init__(feature_definition, *spheres)

        self.__score = 0
        self.__projected_value = 0
        self.__label = 0
        self.__point = None
        self.__projected = None
        self.__projected_identifier = None
        if crystal:
            self.crystal = crystal
            self.coordinate_atm_dict = {_coordinate_str(a.coordinates): a for a in self.crystal.molecule.atoms}
            self.__point = self.get_point()
            self.__projected = self.get_projected()

    @property
    def score(self):
        return self.__score

    @score.setter
    def score(self, value):
        self.__score = value

    @property
    def projected_value(self):
        return self.__projected_value

    @projected_value.setter
    def projected_value(self, value):
        self.__projected_value = value

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, value):
        self.__label = value

    @property
    def point(self):
        return self.__point

    @point.setter
    def point(self, sphere):
        self.__point = [sphere]

    @property
    def projected(self):
        return self.__projected

    @projected.setter
    def projected(self, sphere):
        self.__projected = [sphere]

    @property
    def projected_identifier(self):
        return self.__projected_identifier

    @projected_identifier.setter
    def projected_identifier(self, ident):
        self.__projected_identifier = ident

    # CLASS METHODS

    def get_point(self):
        s = [sphere for sphere in self.spheres if _coordinate_str(sphere.centre) in self.coordinate_atm_dict]

        if len(s) == 0:
            return [self.spheres[0]]
        else:
            return s

    def get_projected(self):
        if len(self.spheres) == 1:
            return None
        else:
            s = [sphere for sphere in self.spheres if _coordinate_str(sphere.centre) in self.coordinate_atm_dict]

            if len(s) == 0:
                return [self.spheres[1]]
            else:
                return [sphere for sphere in self.spheres
                        if not _coordinate_str(sphere.centre) in self.coordinate_atm_dict]

    def to_pymol_str(self, i=0, label=True, transparency=0.8):
        pymol_out = ""
        point_colour = rgb_to_decimal(self.colour)
        point_colour = utilities.Colour(point_colour[0], point_colour[1], point_colour[2], transparency)
        feat_ID = f"{self.identifier}_{i}"
        group = []

        coords = (self.point[0].centre[0],
                  self.point[0].centre[1],
                  self.point[0].centre[2])

        radius = self.point[0].radius

        point_objname = f"{self.identifier}_point_{i}"
        pymol_out += PyMOLCommands.sphere(point_objname, point_colour, coords, radius=radius)
        pymol_out += PyMOLCommands.load_cgo(point_objname, f"{point_objname}_obj")
        group.append(f"{point_objname}_obj")

        if self.score != 0 and label:
            score_name = f"{point_objname}_score"
            pymol_out += PyMOLCommands.pseudoatom(objname=score_name,
                                                  coords=coords,
                                                  label=f'{round(self.score, 1)}')
            group.append(score_name)

        if self.projected:
            proj_coords = (self.projected[0].centre[0],
                           self.projected[0].centre[1],
                           self.projected[0].centre[2])

            projection_objname = f"{self.identifier}_projection_{i}"
            line_objname = f"{self.identifier}_line_{i}"

            group.extend([projection_objname, line_objname])
            if self.projected_value != 0 and label:
                proj_score_name = f"{point_objname}_proj_score"
                pymol_out += PyMOLCommands.pseudoatom(objname=proj_score_name,
                                                      coords=proj_coords,
                                                      label=f'{round(self.projected_value, 1)}')
                group.append(proj_score_name)
            projected_colour = adjust_lightness(self.colour, percentage=30)
            projected_colour = utilities.Colour(projected_colour[0],
                                                projected_colour[1],
                                                projected_colour[2],
                                                transparency)
            projected_radius = self.projected[0].radius

            pymol_out += PyMOLCommands.sphere(projection_objname, projected_colour, proj_coords,
                                              radius=projected_radius)
            pymol_out += PyMOLCommands.load_cgo(projection_objname, f"{projection_objname}_obj")
            pymol_out += PyMOLCommands.line(line_objname, coords, proj_coords, rgb=projected_colour)

        pymol_out += PyMOLCommands.group(feat_ID, group)
        return pymol_out


Pharmacophore.Feature = Feature

if __name__ == "__main__":
    print("hi")
