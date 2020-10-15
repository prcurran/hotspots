import unittest

from hotspots.hs_io import HotspotReader
from ccdc.utilities import _csd_required, PushDir
from ccdc import io
from ccdc.pharmacophore import Pharmacophore, GeometricDescriptors
from hotspots.protein_extension import Protein
from hotspots.grid_extension import Grid
from hotspots.wrapper_pymol import PyMOLFile, PyMOLCommands
from hotspots.pharmacophore_extension \
    import LigandPharmacophoreModel,  InteractionPharmacophoreModel, \
    PharmacophoreModel, ProteinPharmacophoreModel, HotspotPharmacophoreModel, Feature, \
    create_consensus, _create_grids, _features_to_grid, closest_peak_index


import os

_csd = None


def csd():
    global _csd
    if _csd is None:
        _csd = io.EntryReader('csd')
    return _csd

##################################################################################################################

class NoClassMethods(unittest.TestCase):
    def setUp(self):
        self.parent_dir = "testdata/pharmacophore_extension/PharmacophoreModel"
        self.fnames = ["1dmt_ligand.cm", "1r1h_ligand.cm", "1r1j_ligand.cm", "1y8j_ligand.cm"]
        self.pharmacophores = [PharmacophoreModel.from_file(os.path.join(self.parent_dir, f)) for f in self.fnames]

        self.cm_dir = os.path.dirname(os.path.dirname(io.csd_directory()))
        Pharmacophore.read_feature_definitions(os.path.join(self.cm_dir, "CSD_CrossMiner/feature_definitions"))

    def test_create_grids(self):
        grid_dic = _create_grids(self.pharmacophores)
        # Point is unhashable, set comprehension is not possible
        features = [(s.centre[0], s.centre[1], s.centre[2])
                    for p in self.pharmacophores for f in p.features for s in f.spheres]
        # Point causes problems with Grid.value_at_point()
        self.assertTrue(all([grid_dic["ring_planar_projected"].contains_point(p) for p in features]))

    def test_features_to_grid(self):
        pm = PharmacophoreModel()
        pm.feature_definitions = ["acceptor"]
        pts = [[1.0, 1.0, 1.0],
               [1.0, 1.0, 1.0],
               [2.5, 2.5, 2.5],
               [3.5, 3.5, 3.5]]

        g = Grid.initalise_grid(pts, padding=3)
        features = [Feature(pm.feature_definitions["acceptor"],
                            GeometricDescriptors.Sphere(centre=p, radius=1)) for p in pts]

        for f in features:
            f.point = f.spheres[0]

        h = _features_to_grid(features, g)

        # should be 3 peaks, all with a value of 1.0
        self.assertEqual(3, len(h.get_peaks(min_distance=1, cutoff=0)))
        self.assertTrue(all([h.value_at_point(peak) == 1.0 for peak in h.get_peaks(min_distance=1, cutoff=0)]))

    # def test_closest_peak(self):
    #     x = 1

    def test_consensus(self):
        feats = create_consensus(self.pharmacophores)
        self.assertTrue(2, len(feats))

##################################################################################################################

class TestProteinPharmacophoreModel(unittest.TestCase):
    def setUp(self):
        top = os.getcwd()
        # self.prot = Protein.from_file("testdata/6y2g_A/protein.pdb")
        # self.mol_grid = Grid.from_file("testdata/6y2g_A/molA.grd")

        self.binding_site = Protein.from_file("testdata/6y2g_A/binding_site.pdb")
        self.protein_pharmacophore = ProteinPharmacophoreModel()

        self.protein_pharmacophore.feature_definitions = ["acceptor_projected",
                                                          "donor_projected",
                                                          "donor_ch_projected"]

    def test_detect_from_prot(self):
        # create binding site
        # bs = Protein.BindingSiteFromGrid(self.prot, self.mol_grid)
        # self.binding_site = self.prot.copy()
        # remove = {r.identifier for r in self.prot.residues} - {r.identifier for r in bs.residues}
        # for r in remove:
        #     self.binding_site.remove_residue(r)
        # with io.MoleculeWriter("testdata/6y2g_A/binding_site.pdb") as w:
        #     w.write(self.binding_site)
        ###

        self.protein_pharmacophore.detect_from_prot(self.binding_site)

        print(self.protein_pharmacophore.feature_definitions)
        donor_feats = [f for f in self.protein_pharmacophore.detected_features if f.identifier == "donor_projected"]
        self.assertEqual(26, len(donor_feats))
        self.protein_pharmacophore.pymol_visulisation(outdir="testdata/pharmacophore_extension/ProteinPharmacophoreModel/from_prot")


##################################################################################################################

class TestHotspotPharmacophoreModel(unittest.TestCase):
    def setUp(self) -> None:
        x = 1
        with HotspotReader("testdata/pharmacophore_extension/provided_data/out.zip") as r:
            self.hr = [hr for hr in r.read() if hr.identifier == "best_volume"][0]

        # smoothing is really important to this workflow
        for p, g in self.hr.super_grids.items():
            h = g.max_value_of_neighbours()
            h = h.gaussian()
            self.hr.super_grids[p] = h

    def test_noprojections(self):
        p = HotspotPharmacophoreModel()
        p.from_hotspot(self.hr, projections=False)

        top_feats = p.top_features(3)
        p.detected_features = top_feats

        p.pymol_visulisation("testdata/pharmacophore_extension/HotspotPharmacophoreModel/no_projections")

    def test_projections(self):
        p = HotspotPharmacophoreModel()
        p.from_hotspot(self.hr, projections=True)

        p.pymol_visulisation("testdata/pharmacophore_extension/HotspotPharmacophoreModel/projections")



##################################################################################################################
@_csd_required
class TestLigandPharmacophoreModel(unittest.TestCase):
    def setUp(self):
        self.csd = csd()
        self.crystal = csd().crystal('AABHTZ')
        self.crystal.molecule.add_hydrogens()
        self.ligand_pharmacophore = LigandPharmacophoreModel()

    def testSetters(self):
        self.assertEqual(0, len(self.ligand_pharmacophore.features))
        self.ligand_pharmacophore.feature_definitions = ["acceptor"]
        self.assertEqual(1, len(self.ligand_pharmacophore.feature_definitions))

    def testdetect_from_ligand(self):
        wrkdir = "testdata/pharmacophore_extension/LigandPharmacophoreModel/from_ligand"
        with PushDir(wrkdir):
            self.ligand_pharmacophore.feature_definitions = ["acceptor"]
            print(self.ligand_pharmacophore.feature_definitions["acceptor"].point_generator_names)
            self.ligand_pharmacophore.detect_from_ligand(ligand=self.crystal)
            self.assertEqual(5, len(self.ligand_pharmacophore.detected_features))
            # test score setter
            self.assertEqual(0, self.ligand_pharmacophore.detected_features[0].score)
            self.ligand_pharmacophore.detected_features[0].score = 5
            self.assertEqual(5, self.ligand_pharmacophore.detected_features[0].score)

            # test write
            for f in self.ligand_pharmacophore.detected_features:
                self.ligand_pharmacophore.add_feature(f)

            self.ligand_pharmacophore.write(pharmacophore_path="pharmacophore.cm")

            read_me = LigandPharmacophoreModel.from_file("pharmacophore.cm")

            self.assertEqual(len(self.ligand_pharmacophore.features), len(read_me.features))

            print(read_me.features[0].spheres[0].centre)

        # self.ligand_pharmacophore.pymol_visulisation()

    def testto_pymol_str(self):
        self.ligand_pharmacophore.feature_definitions = ["acceptor_projected"]
        self.ligand_pharmacophore.detect_from_ligand(ligand=self.crystal)

        f = PyMOLFile()
        f.commands += self.ligand_pharmacophore.detected_features[0].to_pymol_str()
        f.write("testdata/pharmacophore_extension/LigandPharmacophoreModel/from_ligand/feature_write.py")

    def testdetect_from_pdb(self):
        testpdb = "2vta"
        testhetid = "LZ1"
        testchainid = "A"
        self.ligand_pharmacophore.feature_definitions = ["donor_projected",
                                                         "acceptor_projected"]
        self.ligand_pharmacophore.detect_from_pdb(pdb=testpdb,
                                                  hetid=testhetid,
                                                  chainid=testchainid)

        self.assertEqual(2, len(self.ligand_pharmacophore.detected_features))
        # self.ligand_pharmacophore.pymol_visulisation("testdata/pharmacophore_extension/LigandPharmacophoreModel/from_pdb")

    def testdetect_from_ligand_ensemble(self):
        wrk_dir = "testdata/pharmacophore_extension/LigandPharmacophoreModel/from_ligand_ensemble"
        with PushDir(wrk_dir):
            test_overlay = io.MoleculeReader("test_overlay.mol2")
            ligand_pharmacophore = LigandPharmacophoreModel()
            ligand_pharmacophore.feature_definitions = ["ring_planar_projected"]

            ligand_pharmacophore.detect_from_ligand_ensemble(ligands=test_overlay, cutoff=2)
            # ligand_pharmacophore.pymol_visulisation(outdir="")

            self.assertEqual(2, len(ligand_pharmacophore.detected_features))

    def testdetect_from_ligand_ensemble_cdk2(self):
        wrk_dir = "testdata/pharmacophore_extension/LigandPharmacophoreModel/from_ligand_ensemble_big_all"
        with PushDir(wrk_dir):
            test_overlay = io.MoleculeReader("cdk2_ligands.mol2")
            ligand_pharmacophore = LigandPharmacophoreModel()
            ligand_pharmacophore.feature_definitions = ["ring_planar_projected",
                                                        "donor_projected",
                                                        "acceptor_projected"]

            ligand_pharmacophore.detect_from_ligand_ensemble(ligands=test_overlay, cutoff=2)

            feature_count = 4
            selected = ligand_pharmacophore.top_features(num=feature_count)
            ligand_pharmacophore.detected_features = selected

            self.assertEqual(feature_count, len(ligand_pharmacophore))
            # ligand_pharmacophore.pymol_visulisation(outdir="")

##################################################################################################################

#
#
# class TestInteractionPharmacophoreModel2vta(unittest.TestCase):
#     def setUp(self):
#         self.protein_path = "testdata/pharmacophore_extension/2vta/test_2vta_protonated.pdb"
#         self.hetid = "LZ1"
#         self.chain = "A"
#         self.prot_lig_pharmacophore = InteractionPharmacophoreModel()
#         self.featdefs = ["donor_projected",
#                        "donor_ch_projected",
#                        "acceptor_projected",
#                        "ring"]
#
#         self.prot_lig_pharmacophore.feature_definitions = self.featdefs
#
#     def testligand_pharmacophore(self):
#         self.assertEqual(0, len(self.prot_lig_pharmacophore.features))
#         self.assertEqual(len(self.featdefs), len(self.prot_lig_pharmacophore.feature_definitions))
#
#     def testdetection(self):
#         self.prot_lig_pharmacophore.detect_from_arpeggio(self.protein_path, self.hetid, self.chain)
#         self.assertEqual(6, len(self.prot_lig_pharmacophore.detected_features))
#         self.prot_lig_pharmacophore.pymol_visulisation("testdata/pharmacophore_extension/2vta")
#
#
# class TestInteractionPharmacophoreModel1xkk(unittest.TestCase):
#     def setUp(self):
#         self.protein_path = "testdata/pharmacophore_extension/1xkk/1xkk.pdb"
#         self.hetid = "FMM"
#         self.chain = "A"
#         self.prot_lig_pharmacophore = InteractionPharmacophoreModel()
#         self.prot_lig_pharmacophore.feature_definitions = ["donor_projected",
#                                                          "donor_ch_projected",
#                                                          "acceptor_projected",
#                                                          "ring"]
#
#     # def testligand_pharmacophore(self):
#     #     self.assertEqual(0, len(self.prot_lig_pharmacophore.features))
#     #     self.assertEqual(6, len(self.prot_lig_pharmacophore.feature_definitions))
#
#     def testdetection(self):
#         self.prot_lig_pharmacophore.detect_from_arpeggio(self.protein_path, self.hetid, self.chain)
#         # self.assertEqual(7, len(self.prot_lig_pharmacophore.selected_features))
#         self.prot_lig_pharmacophore.pymol_visulisation("testdata/pharmacophore_extension/1xkk")
#


class TestInteractionPharmacophoreModel1aq1(unittest.TestCase):
    def setUp(self):
        self.protein_path = "/home/pcurran/github_packages/pharmacophores/testdata/alignment/1AQ1_aligned.pdb"
        self.hetid = "STU"
        self.chain = "A"
        self.prot_lig_pharmacophore = InteractionPharmacophoreModel()
        self.prot_lig_pharmacophore.feature_definitions = ["donor_projected",
                                                         "donor_ch_projected",
                                                         "acceptor_projected",
                                                         "ring_planar_projected"]

    # def testligand_pharmacophore(self):
    #     self.assertEqual(0, len(self.prot_lig_pharmacophore.features))
    #     self.assertEqual(6, len(self.prot_lig_pharmacophore.feature_definitions))

    def testdetection(self):
        self.prot_lig_pharmacophore.detect_from_arpeggio(self.protein_path, self.hetid, self.chain)
        # self.assertEqual(7, len(self.prot_lig_pharmacophore.selected_features))
        self.prot_lig_pharmacophore.pymol_visulisation("/home/pcurran/github_packages/pharmacophores/testdata/alignment")



class TestInteractionPharmacophoreModel(unittest.TestCase):
    def setUp(self):
        self.wrkdir = "/home/pcurran/github_packages/pharmacophores/testdata/alignment"

    # def testdetect_interactions(self):
    #     with PushDir(self.wrkdir):
    #         bs = Protein.from_file("1AQ1_bs.pdb")
    #         lig = io.MoleculeReader("1AQ1_STU_aligned.mol2")[0]
    #         ipm = InteractionPharmacophoreModel()
    #
    #         ipm.detect_interactions(bs, lig)

    def testdetect_from_pl_ensemble(self):
        wrkdir = "/home/pcurran/github_packages/pharmacophores/testdata/alignment"
        with PushDir(wrkdir):
            paths = ["1AQ1_aligned.pdb", "1B38_aligned.pdb", "1B39_aligned.pdb", "1CKP_aligned.pdb"]
            hetids = ["STU", "ATP", "ATP", "PVB"]
            chains = ["A", "A", "A", "A"]

            for i in range(len(paths)):
                ipm = InteractionPharmacophoreModel()
                ipm.feature_definitions = ["donor_projected",
                                             "donor_ch_projected",
                                             "acceptor_projected",
                                             "ring_planar_projected"]

                ipm.detect_from_arpeggio(paths[i], hetids[i], chains[i])

                for feat in ipm.detected_features:
                    ipm.add_feature(feat)

                ipm.write(f"{paths[i].split('_')[0]}_{hetids[i]}.cm")
##################################################################################################################

if __name__ == '__main__':
    unittest.main()
