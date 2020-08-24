import unittest

from ccdc.utilities import _csd_required
from ccdc import io
from hotspots.protein_extension import Protein
from hotspots.grid_extension import Grid
from hotspots.pharmacophore_extension \
    import LigandPharmacophoreModel,  InteractionPharmacophoreModel, PharmacophoreModel, ProteinPharmacophoreModel

import os

_csd = None


def csd():
    global _csd
    if _csd is None:
        _csd = io.EntryReader('csd')
    return _csd


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
        donor_feats = [f for f in self.protein_pharmacophore.selected_features if f.identifier == "donor_projected"]
        self.assertEqual(26, len(donor_feats))
        self.protein_pharmacophore.pymol_visulisation(outdir="testdata/pharmacophore_extension/ProteinPharmacophoreModel/from_prot")


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
        self.ligand_pharmacophore.feature_definitions = ["acceptor"]
        self.ligand_pharmacophore.detect_from_ligand(ligand=self.crystal)
        self.assertEqual(5, len(self.ligand_pharmacophore.selected_features))
        # test score setter
        self.assertEqual(0, self.ligand_pharmacophore.selected_features[0].score)
        self.ligand_pharmacophore.selected_features[0].score = 5
        self.assertEqual(5, self.ligand_pharmacophore.selected_features[0].score)
        self.ligand_pharmacophore.pymol_visulisation("testdata/pharmacophore_extension/LigandPharmacophoreModel/from_ligand")


    def testdetect_from_pdb(self):
        testpdb = "2vta"
        testhetid = "LZ1"
        testchainid = "A"
        self.ligand_pharmacophore.feature_definitions = ["donor_projected", "acceptor_projected"]
        self.ligand_pharmacophore.detect_from_pdb(pdb=testpdb,
                                                  hetid=testhetid,
                                                  chainid=testchainid)
        self.assertEqual(2, len(self.ligand_pharmacophore.features))

    # TODO: Test Methods used in detect_from_ligand_ensemble
    def testdetect_from_ligand_ensemble(self):
        test_overlay = io.MoleculeReader("testdata/pharmacophore_extension/test_overlay.mol2")
        # # projected feature
        # ligand_pharmacophore_1 = LigandPharmacophoreModel()
        # ligand_pharmacophore_1.feature_definitions = ["donor_projected"]
        # ligand_pharmacophore_1.detect_from_ligand_ensemble(ligands=test_overlay)
        # self.assertEqual(3, len(ligand_pharmacophore_1.selected_features))
        #
        # # non projected feature
        # ligand_pharmacophore_2 = LigandPharmacophoreModel()
        # ligand_pharmacophore_2.feature_definitions = ["hydrophobe"]
        # ligand_pharmacophore_2.detect_from_ligand_ensemble(ligands=test_overlay)
        # self.assertEqual(9, len(ligand_pharmacophore_2.selected_features))

    # Test Methods used in detect_from_pdb_ensemble
        ligand_pharmacophore_3 = LigandPharmacophoreModel()
        ligand_pharmacophore_3.feature_definitions = ["acceptor_projected",
                                                      "donor_projected",
                                                      "ring_planar_projected"]

        ligand_pharmacophore_3.detect_from_ligand_ensemble(ligands=test_overlay, density=0.6)

        ligand_pharmacophore_3.pymol_visulisation(outdir="testdata/pharmacophore_extension")


class TestInteractionPharmacophoreModel2vta(unittest.TestCase):
    def setUp(self):
        self.protein_path = "testdata/pharmacophore_extension/2vta/test_2vta_protonated.pdb"
        self.hetid = "LZ1"
        self.chain = "A"
        self.prot_lig_pharmacophore = InteractionPharmacophoreModel()
        self.prot_lig_pharmacophore.feature_definitions = ["donor_projected",
                                                           "donor_ch_projected",
                                                           "acceptor_projected",
                                                           "ring"]

    def testligand_pharmacophore(self):
        self.assertEqual(0, len(self.prot_lig_pharmacophore.features))
        self.assertEqual(6, len(self.prot_lig_pharmacophore.feature_definitions))

    def testdetection(self):
        self.prot_lig_pharmacophore.detect_from_protein_ligand_complex(self.protein_path, self.hetid, self.chain)
        self.assertEqual(7, len(self.prot_lig_pharmacophore.selected_features))
        self.prot_lig_pharmacophore.pymol_visulisation("testdata/pharmacophore_extension/2vta")


class TestInteractionPharmacophoreModel1xkk(unittest.TestCase):
    def setUp(self):
        self.protein_path = "testdata/pharmacophore_extension/1xkk/1xkk.pdb"
        self.hetid = "FMM"
        self.chain = "A"
        self.prot_lig_pharmacophore = InteractionPharmacophoreModel()
        self.prot_lig_pharmacophore.feature_definitions = ["donor_projected",
                                                         "donor_ch_projected",
                                                         "acceptor_projected",
                                                         "ring"]

    # def testligand_pharmacophore(self):
    #     self.assertEqual(0, len(self.prot_lig_pharmacophore.features))
    #     self.assertEqual(6, len(self.prot_lig_pharmacophore.feature_definitions))

    def testdetection(self):
        self.prot_lig_pharmacophore.detect_from_protein_ligand_complex(self.protein_path, self.hetid, self.chain)
        # self.assertEqual(7, len(self.prot_lig_pharmacophore.selected_features))
        self.prot_lig_pharmacophore.pymol_visulisation("testdata/pharmacophore_extension/1xkk")


class TestA(unittest.TestCase):
    def testsmarts(self):

        p = PharmacophoreModel()
        p.feature_definitions=["donor_projected"]
        dp = p.feature_definitions["donor_projected"]
        for s in dp.SMARTS_definitions:
            print(str(s))


if __name__ == '__main__':
    unittest.main()
