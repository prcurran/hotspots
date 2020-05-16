import unittest
import os

from hotspots.wrapper_arpeggio import Arpeggio
from hotspots.wrapper_arpeggio import _get_protein, _get_ligand
from ccdc.protein import Protein
from ccdc.molecule import Molecule, Atom, Coordinates
from ccdc.io import MoleculeReader, MoleculeWriter

from scipy.spatial import distance
import numpy as np
import shutil


class TestOthers(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = "testdata/wrapper_arpeggio/prepare"
        self.copydir = "testdata/wrapper_arpeggio/prepare/copydir"

        self.examples = {"2vta": "LZ1",     # single chain, single ligand
                         "3img": "BZ2",     # multiple chains, multiple ligands
                         "1xkk": "FMM"}     #

    def test_get_protein(self):
        # 2VTA example
        for pdb, hetid in self.examples.items():
            # from PDB code only
            tmp1 = os.path.join(self.tmpdir, "pdb")
            protein1 = _get_protein(pdb, tmp1)
            self.assertTrue(isinstance(protein1, Protein))                          # GOT A PROTEIN
            self.assertTrue("H" not in {atm.atomic_symbol for atm in protein1.atoms})       # NO H's

            tmp2 = os.path.join(self.tmpdir, "file")
            protein2 = _get_protein(pdb, tmp2, f"testdata/wrapper_arpeggio/prepare/copydir/{pdb}.pdb")
            self.assertTrue(isinstance(protein2, Protein))
            self.assertTrue("H" in {atm.atomic_symbol for atm in protein2.atoms})

    def test_get_ligand(self):
        sources = ["file", "pdb"]
        for source in sources:
            tmp = f"testdata/wrapper_arpeggio/prepare/{source}"
            for pdb, hetid in self.examples.items():
                protein1 = Protein.from_file(os.path.join(tmp, f"{pdb}_clean.pdb"))
                ligand1 = _get_ligand(protein1, hetid, "A")
                self.assertTrue(isinstance(ligand1, Molecule))


class TestArpeggio(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = os.path.abspath("testdata/wrapper_arpeggio/tmpdir")
        source = "file"
        self.data = f"testdata/wrapper_arpeggio/prepare/{source}"

        self.arpeggio_2vta = Arpeggio(pdb_code="2vta",
                                      hetid="LZ1",
                                      chain="A",
                                      tmpdir=os.path.join(self.tmpdir, "lz1"))
        self.testlz1 = MoleculeReader(os.path.join(self.data, "LZ1.mol2"))[0]
        self.lz1_resid = "A/1301/"

        self.arpeggio_1xkk = Arpeggio(pdb_code="1xkk",
                                      hetid="FMM",
                                      chain="A",
                                      tmpdir=os.path.join(self.tmpdir, "fmm"))
        self.testfmm = MoleculeReader(os.path.join(self.data, "FMM.mol2"))[0]
        self.fmm_resid = "A/91/"

        # shutil.rmtree(self.arpeggio_2vta.tmpdir)
        # shutil.rmtree(self.arpeggio_1xkk.tmpdir)

    def testinit(self):
        self.assertEqual(1301, self.arpeggio_2vta.ligand_index)
        self.assertEqual(91, self.arpeggio_1xkk.ligand_index)

    def testrun(self):
        self.arpeggio_2vta.run()
        self.assertTrue(len(os.listdir(self.arpeggio_2vta.tmpdir)) > 1)

        self.arpeggio_1xkk.run()
        self.assertTrue(len(os.listdir(self.arpeggio_1xkk.tmpdir)) > 1)

    def testread_2vta(self):
        extensions = ["sift", "contacts", "atomtypes", "rings", "ri", "ari"]
        fs = set([f for f in os.listdir(self.arpeggio_2vta.tmpdir) if f.split(".")[1] in extensions])

        self.assertEqual(len(extensions), len(fs))

        summary = self.arpeggio_2vta._get_summary()
        self.assertEqual(9, len(summary))
        self.assertEqual(2, len(summary.loc[summary.hbond == 1].hbond.values))     # there are 2 hbond interactions
        self.assertEqual(1, len(summary.loc[summary.aromatic == 1].hbond.values))   # there is 1 aromatic interactions

        contacts = self.arpeggio_2vta._get_contacts()
        self.assertEqual(115, len(contacts))

        atomtypes = self.arpeggio_2vta._get_atom_types()
        self.assertEqual(len(self.arpeggio_2vta.protein.atoms), len(atomtypes))
        # N is aromatic and donor
        self.assertEqual(2, len(atomtypes.loc[atomtypes.atom == f"{self.lz1_resid}N"].atomtype.values[0]))

        rings = self.arpeggio_2vta._get_rings()
        # Indazole is a bicycle ;)
        ligand_rings = rings.loc[rings.resid == self.lz1_resid].ringid.values
        self.assertEqual(2, len(ligand_rings))

        ring_to_ring = self.arpeggio_2vta._get_ring_to_ring()
        # 2 ring to ring contact
        self.assertEqual(2, len(ring_to_ring.loc[(ring_to_ring.resid1 == self.lz1_resid) |
                                                 (ring_to_ring.resid2 == self.lz1_resid)]))

        atom_to_ring = self.arpeggio_2vta._get_atom_to_ring()
        # ALA,LEU above and below 5-ring + another one
        # NB ringid is the same on each run
        self.assertEqual(3, len(atom_to_ring.loc[atom_to_ring.resid == self.lz1_resid]))

    def testread_1xkk(self):
        extensions = ["sift", "contacts", "atomtypes", "rings", "ri", "ari"]
        fs = set([f for f in os.listdir(self.arpeggio_1xkk.tmpdir) if f.split(".")[1] in extensions])

        self.assertEqual(len(extensions), len(fs))

        summary = self.arpeggio_1xkk._get_summary()

        summary.to_csv(os.path.join(self.tmpdir, "test.csv"))
        self.assertEqual(40, len(summary))
        self.assertEqual(5, len(summary.loc[summary.hbond == 1].hbond.values))     # there are 4 hbond interactions
        self.assertEqual(2, len(summary.loc[summary.aromatic == 1].hbond.values))   # there is 2 aromatic interactions

        contacts = self.arpeggio_1xkk._get_contacts()
        self.assertEqual(482, len(contacts))
        contacts.to_csv(os.path.join(self.tmpdir, "test_contacts.csv"))
        atomtypes = self.arpeggio_1xkk._get_atom_types()
        self.assertEqual(len(self.arpeggio_1xkk.protein.atoms), len(atomtypes))

        rings = self.arpeggio_1xkk._get_rings()
        # FMM has 5 rings
        ligand_rings = rings.loc[rings.resid == self.fmm_resid].ringid.values
        self.assertEqual(5, len(ligand_rings))

    def testget_hbonds_2vta(self):
        hbonds = []
        for atm in self.testlz1.heavy_atoms:
            hbonds.extend(self.arpeggio_2vta.get_hbonds(atm, "INTER"))

        self.assertEqual(2, len(hbonds))

        # check the Interaction construction
        hbond = [h for h in hbonds if h.point_identifier == "LZ1/N"][0]

        self.assertTrue(isinstance(hbond.point, Coordinates))
        self.assertTrue(isinstance(hbond.point_atom, Atom))
        self.assertEqual("LZ1/N", hbond.point_identifier)
        self.assertEqual("donor", hbond.point_type)

        self.assertTrue(isinstance(hbond.projected, Coordinates))
        self.assertTrue(isinstance(hbond.point_atom, Atom))
        self.assertEqual("A/81/O", hbond.projected_identifier)
        self.assertEqual("hbond", hbond.interation_type)

    def testget_hbonds_1xkk(self):
        hbonds = []
        for atm in self.testfmm.heavy_atoms:
            hbonds.extend(self.arpeggio_1xkk.get_hbonds(atm))

        self.assertEqual(5, len(hbonds))

        # check the Interaction construction
        hbond = [h for h in hbonds if h.point_identifier == "FMM/O3"][0]
        self.assertTrue(isinstance(hbond.point, Coordinates))
        self.assertTrue(isinstance(hbond.point_atom, Atom))
        self.assertEqual("FMM/O3", hbond.point_identifier)
        self.assertEqual("acceptor", hbond.point_type)

        self.assertTrue(isinstance(hbond.projected, Coordinates))
        self.assertTrue(isinstance(hbond.point_atom, Atom))
        self.assertEqual("A/81/O4", hbond.projected_identifier)
        self.assertEqual("hbond", hbond.interation_type)

    def testget_weak_hbonds_2vta(self):
        weak_hbonds = []
        for atm in self.testlz1.heavy_atoms:
           weak_hbonds.extend(self.arpeggio_2vta.get_weak_hbonds(atm, "INTER"))
        self.assertEqual(2, len(weak_hbonds))

    def testget_weak_hbonds_1xkk(self):
        weak_hbonds = []
        for atm in self.testfmm.heavy_atoms:
           weak_hbonds.extend(self.arpeggio_1xkk.get_weak_hbonds(atm))

        self.assertEqual(15, len(weak_hbonds))

    def testget_ligand_atom_to_ring_bonds(self):
        # self.arpeggio_1xkk = Arpeggio(pdb_code="1xkk",
        #                               hetid="FMM",
        #                               chain="A",
        #                               tmpdir=os.path.join(self.tmpdir, "fmm"))
        #
        # self.arpeggio_1xkk.run()
        #
        # self.testfmm = MoleculeReader("testdata/wrapper_arpeggio/FMM.mol2")[0]
        # TODO: Find an example with some ligand atom to protein ring contacts
        aromatic_contacts = []
        for atm in self.testlz1.heavy_atoms:
            aromatic_contacts.extend(self.arpeggio_2vta.get_ligand_atom_to_ring_bonds(atm))

        self.assertEqual(0, len(aromatic_contacts))

    def testget_ligand_ring_to_atom_bonds_2vta(self):
        aromatic_contacts = []
        rings = self.arpeggio_2vta._get_rings()
        ligand_rings = rings.loc[rings.resid == self.lz1_resid].ringid.values
        for r in ligand_rings:
            aromatic_contacts.extend(self.arpeggio_2vta.get_ligand_ring_to_atom_bonds(r))

        self.assertEqual(3, len(aromatic_contacts))

    def testget_ligand_ring_to_atom_bonds_1xkk(self):
        aromatic_contacts = []
        rings = self.arpeggio_1xkk._get_rings()
        ligand_rings = rings.loc[rings.resid == self.fmm_resid].ringid.values
        for r in ligand_rings:
            aromatic_contacts.extend(self.arpeggio_1xkk.get_ligand_ring_to_atom_bonds(r))

        self.assertEqual(11, len(aromatic_contacts))

    def testget_ring_to_ring_bonds_2vta(self):
        aromatic_contacts = []
        rings = self.arpeggio_2vta._get_rings()
        ligand_rings = rings.loc[rings.resid == self.lz1_resid].ringid.values
        for r in ligand_rings:
            aromatic_contacts.extend(self.arpeggio_2vta.get_ring_to_ring_bonds(r))

        self.assertEqual(1, len(aromatic_contacts))

        aro = aromatic_contacts[0]
        self.assertTrue(isinstance(aro.point, Coordinates))
        self.assertTrue(aro.point_atom is None)
        self.assertEqual(f"LZ1/ringid_{ligand_rings[0]}", aro.point_identifier)

        self.assertTrue(isinstance(aro.projected, Coordinates))
        self.assertTrue(aro.projected_atom is None)
        # self.assertEqual("?", aro.projected_identifier)
        self.assertEqual("aromatic", aro.interation_type)

    def testget_ring_to_ring_bonds_1xkk(self):
        aromatic_contacts = []
        rings = self.arpeggio_2vta._get_rings()
        ligand_rings = rings.loc[rings.resid == self.fmm_resid].ringid.values
        for r in ligand_rings:
            aromatic_contacts.extend(self.arpeggio_1xkk.get_ring_to_ring_bonds(r))

        self.assertEqual(1, len(aromatic_contacts))

        aro = aromatic_contacts[0]
        self.assertTrue(isinstance(aro.point, Coordinates))
        self.assertTrue(aro.point_atom is None)
        self.assertEqual(f"FMM/ringid_{ligand_rings[0]}", aro.point_identifier)

        self.assertTrue(isinstance(aro.projected, Coordinates))
        self.assertTrue(aro.projected_atom is None)
        # self.assertEqual("?", aro.projected_identifier)
        self.assertEqual("aromatic", aro.interation_type)

    def testcreate_feature_list(self):
        atom_features, ring_features = self.arpeggio_2vta.create_feature_list("INTER")
        fs = atom_features + ring_features
        pairs = [[[f.point.x, f.point.y, f.point.z], [f.projected.x, f.projected.y, f.projected.z]] for f in fs]
        d = [distance.euclidean(point, proj) for point, proj in pairs]
        self.assertTrue(all(np.less(d, 8)))


if __name__ == '__main__':
    unittest.main()
