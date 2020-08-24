import os
import tempfile
import colorsys

from scipy.spatial import distance
import numpy as np

from ccdc.pharmacophore import Pharmacophore
from ccdc.io import csd_directory, CrystalReader, MoleculeWriter
from ccdc.molecule import Molecule, Atom, Coordinates
from ccdc.protein import Protein
from ccdc.descriptors import GeometricDescriptors
from ccdc import utilities

from hotspots.grid_extension import Grid
from hotspots.wrapper_pdb import PDBResult
from hotspots.wrapper_arpeggio import Arpeggio
from hotspots.template_strings import pymol_imports, pymol_sphere, pymol_group, pymol_line, \
    pymol_isosurface, pymol_pseudoatom, pymol_set_color, pymol_load

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
    return np.array([t[0], t[1], t[2]])


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
        self.__ligands = None
        self.__protein = None
        self.__selected_features = None
        self.__feature_point_grids = None

        # self.interaction_partners = {"donor_projected": {"hbond": ["acceptor_projected"],
        #                                                  "aromatic": ["ring", "ring_planar_projected"]
        #                                                  },
        #
        #                              "acceptor_projected": {"hbond":["donor_projected"],
        #                                                     "weak_hbond":["donor_ch_projected"]
        #                                                     },
        #
        #                              }
        #
        # ["donor_projected",
        #  "donor_ch_projected",
        #  ,
        #  "ring_planar_projected",
        #  "ring",
        #  "hydrophobe"]

    def __len__(self):
        return len(self.__features)

    @property
    def ligands(self):
        return self.__ligands

    @ligands.setter
    def ligands(self, ligands):
        self.__ligands = ligands

    @property
    def protein(self): return self.__protein

    @protein.setter
    def protein(self, protein): self.__protein = protein

    @property
    def feature_point_grids(self):
        return self.__feature_point_grids

    @feature_point_grids.setter
    def feature_point_grids(self, dict):
        self.__feature_point_grids = dict

    @property
    def selected_features(self):
        return self.__selected_features

    @selected_features.setter
    def selected_features(self, features):
        self.__selected_features = features

    @property
    def feature_definitions(self): return self.__feature_definitions

    @feature_definitions.setter
    def feature_definitions(self, feature_types):
        # reset before choosing
        self.__feature_definitions = {k: v for k, v in self.__feature_options.items()
                                      if any([k == ft for ft in feature_types])}

    # CLASS METHODS
    def features_to_pymol_strings(self, features):
        """
        creates the code to visualise `ccdc.pharmacophore.Pharmacophore.Features` in PyMOL

        :param list of (`ccdc.pharmacophore.Pharmacophore.Features`) features: features to be visualised
        :return: str python code
        """
        group_dic = {ident: {} for ident in self.feature_definitions.keys()}
        pymol_out = ''
        for i, feat in enumerate(features):
            identifier = feat.identifier
            point_colour = rgb_to_decimal(feat.colour)
            feat_ID = f"{identifier}_{i}"
            coords = [feat.point[0].centre[0],
                      feat.point[0].centre[1],
                      feat.point[0].centre[2]]

            point_objname = f"{identifier}_point_{i}"
            group_dic[identifier].update({feat_ID: [point_objname]})
            pymol_out += PyMOLCommands.sphere(point_objname, point_colour, coords, radius=0.5)
            pymol_out += PyMOLCommands.load_cgo(point_objname, f"{point_objname}_obj")

            if feat.projected:
                proj_coords = [feat.projected[0].centre[0],
                               feat.projected[0].centre[1],
                               feat.projected[0].centre[2]]

                projection_objname = f"{identifier}_projection_{i}"
                line_objname = f"{identifier}_line_{i}"
                group_dic[identifier][feat_ID].extend([projection_objname, line_objname])

                projected_colour = adjust_lightness(feat.colour, percentage=30)

                pymol_out += PyMOLCommands.sphere(projection_objname, projected_colour, proj_coords, radius=0.5)
                pymol_out += PyMOLCommands.load_cgo(projection_objname, f"{projection_objname}_obj")
                pymol_out += PyMOLCommands.line(line_objname, coords, proj_coords, rgb=projected_colour)
                if feat.projected_identifier and self.protein:
                    resnum = feat.projected_identifier.split("/")[1]
                    pymol_out += f'\ncmd.select("sele", "resi {resnum}")\ncmd.show("sticks", "sele")'

            pymol_out += PyMOLCommands.group(feat_ID, group_dic[identifier][feat_ID])
            pymol_out += PyMOLCommands.group(f"{identifier}_pts", group_dic[identifier].keys())
            # !!! BE CAREFUL NOT TO REUSE PYMOL NAMESPACES, IT WILL END IN TEARS !!!
        pymol_out += PyMOLCommands.group("ligand_pharmacophore", [f"{a}_pts" for a in group_dic.keys()])

        return pymol_out

    def pymol_visulisation(self, outdir=None):
        if not outdir:
            outdir = os.getcwd()

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        if self.ligands:
            with MoleculeWriter(os.path.join(outdir, "ligands.mol2")) as w:
                for ligand in self.ligands:
                    w.write(ligand.molecule)

        if self.protein:
            with MoleculeWriter(os.path.join(outdir, "protein.mol2")) as w:
                w.write(self.protein.molecule)

        self.pymol_out = PyMOLFile()
        # with open(os.path.join(outdir, "pymol_file.py"), "w") as pymol_file:
        #     pymol_out = pymol_imports()

        if self.ligands:
            self.pymol_out.commands += PyMOLCommands.load("ligands.mol2", "ligands")

        if self.protein:
            self.pymol_out.commands += PyMOLCommands.load("protein.mol2", "protein")

        # write out point spheres and projection sphere and lines if applicable
        self.pymol_out.commands += self.features_to_pymol_strings(self.selected_features)

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

            self.pymol_out.commands += PyMOLCommands.group("grds", self.feature_point_grids.keys())
            self.pymol_out.commands += PyMOLCommands.group("grds", [f"surface_{a}" for a in self.feature_point_grids])

        self.pymol_out.write(os.path.join(outdir, "pymol_file.py"))

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
        self.selected_features = []
        for fd in self.feature_definitions.values():
            detected_feats = fd.detect_features(prot)
            if len(detected_feats) != 0:
                for f in detected_feats:
                    self.selected_features.append(f)


class LigandPharmacophoreModel(PharmacophoreModel):
    """
    Derived Class for representing pharmacophores generated from ligand atom types. Atom
    typing is carried out using built-in CrossMiner SMARTS definitions.

        - Pharmacophore generated from PDB codes use the RCSB PDB RESTful web service
        (https://www.rcsb.org/pdb/software/rest.do)

        - Ligand ensembles points are generated by collecting points on a grid and peak picking

    """
    class FeatureCluster:
        def __init__(self, point, features, feature_def, cutoff):
            self.__point = point
            self.__features = features
            self.__feature_def = feature_def

            self._cutoff = cutoff

            # make from projections
            self.__projection_grid = self.create_projection_grid()
            if not self.__projection_grid:
                self.__projection_grid_peaks = []
            else:
                self.__projection_grid_peaks = self.find_projection_peaks()
            self.__new_features = self.create_new_features()

        def __repr__(self):
            return f"Peak: {_coordinate_str(self.__point, 2)}"

        def __len__(self):
            return len(self.__features)

        # PROPERTIES
        @property
        def point(self): return self.__point

        @property
        def features(self): return self.__features

        @property
        def projection_grid(self): return self.__projection_grid

        @property
        def projection_peaks(self): return self.__projection_grid_peaks

        @property
        def feature_def(self): return self.__feature_def

        @property
        def new_features(self): return self.__new_features

        # CLASS METHODS
        def create_projection_grid(self):
            """
            create the projected spheres for the summary pharmacophore features

            :return: `hotspots.grid_extension.Grid`
            """
            # Do the features have projections ?
            if len([f for f in self.features if f.projected]) == 0:
                return None
            else:
                g = Grid.initalise_grid(coords=[Coordinates(x=feature.projected[0].centre[0],
                                                            y=feature.projected[0].centre[1],
                                                            z=feature.projected[0].centre[2])
                                                for feature in self.features],
                                        spacing=0.1,
                                        padding=2)

                for f in self.features:
                    if f.projected:
                        px, py, pz = [f.projected[0].centre[0],
                                      f.projected[0].centre[1],
                                      f.projected[0].centre[2]]
                        g.set_sphere(point=[px, py, pz],
                                     value=1,
                                     radius=2,
                                     scaling='linear')
                return g

        def find_projection_peaks(self):
            self.proj_by_score = {peak: self.projection_grid.value_at_point(peak)
                                  for peak in self.projection_grid.get_peaks(min_distance=4, cutoff=0)}
            accepted = [k for k, v in self.proj_by_score.items() if round(v) >= self._cutoff]
            return accepted

        def create_new_features(self):
            """
            create new features from points and projections over the cutoff value
            :return: list of (`hotspots.pharmacophore_extension.Feature`)
            """
            new_feats = []
            point = GeometricDescriptors.Sphere(centre=self.point, radius=1)
            # prevent double counting of features which have > 1 projection e.g. SP3 donors (NH3+)
            #       NB: self.projection has already been filtered, use self.proj_by_score instead
            if len(self.projection_peaks) == 0:
                score = float(len(self.features))
            else:
                score = round(float(float(len(self.features))/len(self.proj_by_score)), 2)

            if score >= self._cutoff:
                if len(self.projection_peaks) == 0:
                    feat = Pharmacophore.Feature(self.feature_def, point)
                    feat.point = point
                    feat.score = score
                    new_feats.append(feat)

                else:
                    for a in self.projection_peaks:
                        print(a)
                        proj = GeometricDescriptors.Sphere(centre=a, radius=1)
                        print(proj)
                        feat = Pharmacophore.Feature(self.feature_def, point, proj)
                        feat.point = point
                        feat.projected = proj
                        print(feat.projected)
                        feat.score = score
                        new_feats.append(feat)

            return new_feats

        # STATIC METHODS
        @staticmethod
        def cluster(all_features, all_peaks, feature_def, cutoff, max_distance=2):
            """
            For each peak, the closest feature point is assigned as a peak member

            :param tuple peak: coordinates x, y, z
            :param `Feature`: a pharmacophore feature
            :return: None
            """
            feature_by_peak = {}
            peak_id_by_coord = {_coordinate_str(p): p for p in all_peaks}
            for feature in all_features:
                f = np.array([[feature.point[0].centre[0],
                               feature.point[0].centre[1],
                               feature.point[0].centre[2]]])

                peaks_array = np.array(all_peaks)

                d = distance.cdist(f, peaks_array)
                index = np.argmin(d)

                if d[0][index] < max_distance:
                    try:
                        feature_by_peak[_coordinate_str(all_peaks[index])].append(feature)

                    except:
                        feature_by_peak.update({_coordinate_str(all_peaks[index]): [feature]})
            return [LigandPharmacophoreModel.FeatureCluster(point=peak_id_by_coord[k],
                                                            features=v,
                                                            feature_def=feature_def,
                                                            cutoff=cutoff)
                    for k, v in feature_by_peak.items()]

    def __init__(self):
        super().__init__()

    # from PDB
    def _align_proteins(self):
        return 0

    def _cluster_ligands(self):
        return 0

    def add_self_features(self):
        # detection and addition separated to enable the selection to be adjusted after the
        # inital detection
        for feature in self.features:
            print(f"adding {feature.identifier} feature")
            self.add_feature(feature)

    def detect_from_ligand(self, ligand):
        """
        Add features to pharmacophore model from a ligand

        :param `ccdc.crystal.Crystal` or `ccdc.molecule.Molecule` ligand: a ligand
        :return:
        """
        if isinstance(ligand, Molecule):
            ligand = self._get_crystal(ligand)

        self.ligands = [ligand]

        self.selected_features = []
        for fd in self.feature_definitions.values():
            detected_feats = fd.detect_features(ligand)
            if len(detected_feats) != 0:
                for f in detected_feats:
                    self.selected_features.append(f)

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

    def detect_from_ligand_ensemble(self, ligands, density=0.5):
        """
        detect features (included projected points) from a series of overlaid ligands
        :param `ccdc.crystal.Crystal`, `ccdc.molecule.Molecule` ligands: collection of 1 or more ligands
        :return list of (`ccdc.pharmacophore.Pharmacophore.Feature`): list of new features
        """
        new_features = []

        if isinstance(ligands[0], Molecule):
            ligands = [self._get_crystal(ligand) for ligand in ligands]

        self.ligands = ligands

        g = Grid.initalise_grid(coords=[a.coordinates for l in ligands for a in l.molecule.atoms], spacing=0.1)
        self.feature_point_grids = {k: g.copy() for k in self.feature_definitions.keys()}

        for identifier, fd in self.feature_definitions.items():
            # detect features in ligands
            all_features = [feature for ligand in ligands for feature in fd.detect_features(ligand)]
            # add features to a grid
            self.feature_point_grids[identifier] = self._features_to_grid(all_features,
                                                                          self.feature_point_grids[identifier])
            all_peaks = self.feature_point_grids[identifier].get_peaks(min_distance=4, cutoff=0)
            # Feature below cutoff removed
            cutoff = len(self.ligands) * density
            # assign projected spheres to point sphere, remove projections below cutoff.
            feature_clusters = LigandPharmacophoreModel.FeatureCluster.cluster(all_features, all_peaks, fd, cutoff)

            new_features.extend([feat for fc in feature_clusters for feat in fc.new_features])

        self.selected_features = new_features

        return new_features

    def from_pdb_ensemble(self):
        # cut and paste from odl code
        x = 1

    @staticmethod
    def _features_to_grid(all_features, grid):
        seen = []
        for feature in all_features:
            x, y, z = [feature.point[0].centre[0], feature.point[0].centre[1], feature.point[0].centre[2]]
            # PREVENT DOUBLE COUNTING
            # for each feature definition, only add the FeaturePoint (point of corresponding atom) once
            # for example, C=O Oxygen will generate 2 acceptor_projected features, one for each lone pair
            # we only want to count the oxygen atom position once. (It can be unpacked to 2 features later)
            coordstr = _coordinate_str(feature.point[0].centre)
            if not coordstr in seen:
                grid.set_sphere(point=[x, y, z], value=1, radius=1, scaling='linear')
                seen.append(coordstr)
        return grid


class HotspotPharmacophoreModel(PharmacophoreModel):
    """
    Derived Class for representing pharmacophores generated from Fragment Hotspot Maps.
    """
    def __init__(self):
        super().__init__()

    def from_hotspot(self):
        x = 1


class InteractionPharmacophoreModel(PharmacophoreModel):
    """
    Derived Class for representing pharmacophores generated from Arpeggio.

    This approach uses detected protein-ligand interactions to generate pharmacophore features
    """
    def __init__(self):
        super().__init__()
    #
    # def detect_from_pl_ensemble(self, protein_paths, hetids, chains):

    def detect_from_protein_ligand_complex(self, protein_path, hetid, chain):
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
        assert("H" in {atom.atomic_symbol for atom in protein.atoms[:50]})

        # assertion 2: protein must contain the ligand of interest
        assert(len([l for l in protein.ligands if l.identifier.split(":")[0] == chain and
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

                new = Feature(fd, point, projected)
                new.point = point
                new.point_identifier = interaction_features[i].point_identifier
                new.projected = projected
                new.projected_identifier = interaction_features[i].projected_identifier
                print(new.projected_identifier)

                new_features.append(new)

        self.selected_features = new_features
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
    def score(self): return self.__score

    @score.setter
    def score(self, value): self.__score = value

    @property
    def label(self): return self.__label

    @label.setter
    def label(self, value): self.__label = value

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
    def projected_identifier(self): return self.__projected_identifier

    @projected_identifier.setter
    def projected_identifier(self, ident): self.__projected_identifier = ident

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


Pharmacophore.Feature = Feature


if __name__ == "__main__":
    print("hi")