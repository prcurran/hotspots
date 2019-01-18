import os
import argparse
import tempfile

from ccdc.protein import Protein
from ccdc.descriptors import MolecularDescriptors
from ccdc.io import MoleculeReader

from pdb_python_api import PDBResult


class Organiser(argparse.ArgumentParser):
    """
    class to handle docking run with hotspot insights
    """

    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        # handle command line arguments
        self.add_argument(
            'protein',
            help='pdb_code of protein which was used in docking'
        )

        self.add_argument(
            'reference',
            help='pdb_code of reference'
        )

        self.add_argument(
            'chemical_id',
            help='PDB identifier for the docked ligand'
        )

        self.add_argument(
            'results',
            help='path to results files'
        )

        self.add_argument(
            '-r', '--chain_ref',
            default='A',
            help='Chain to used for alignment'
        )

        self.add_argument(
            '-p', '--chain_protein',
            default='A',
            help='Chain to used for alignment'
        )
        self.args = self.parse_args()
        self.tmp = tempfile.mkdtemp()

        # download protein
        PDBResult(self.args.protein).download(self.tmp)
        self.protein = Protein.from_file(os.path.join(self.tmp,
                                                      self.args.protein + ".pdb"))
        self.protein.add_hydrogens()

        # download reference
        PDBResult(self.args.reference).download(self.tmp)
        ref = Protein.from_file(os.path.join(self.tmp,
                                             self.args.reference + ".pdb"))
        ref.add_hydrogens()

        self.ref = self._align(self.protein, ref)
        self.reference_ligand = self._extract_ligands(protein=self.ref,
                                                      ligand=self.args.chemical_id,
                                                      chain=self.args.chain_ref)[0]

        self.results = MoleculeReader(os.path.join(os.path.dirname(os.path.realpath(__file__)), self.args.results))

    def rmsd(self):
        """
        prepares molecules for rmsd calculation
        :param a:
        :param b:
        :return:
        """
        self.reference_ligand = self._clean_mol(self.reference_ligand)
        return[MolecularDescriptors.rmsd(self.reference_ligand, self._clean_mol(l)) for l in self.results]

    @staticmethod
    def _clean_mol(mol):
        """
        cleans molecule
        removes metal centers, hydrogens and dummy atoms (added during docking)
        :return:
        """
        for atm in mol.atoms:
            if atm.is_metal:
                mol.remove_atom(atm)
        mol.remove_hydrogens()
        mol.remove_unknown_atoms()
        return mol

    @staticmethod
    def _extract_ligands(protein, ligand, chain):
        """
        remove ligands from protein
        :return:
        """
        ligs = []
        for lig in protein.ligands:
            if lig.identifier.split(":")[0] == chain and lig.identifier.split(":")[1][0:3] == ligand:
                ligs.append(lig)
        return ligs

    @staticmethod
    def _align(reference, target, reference_chain='A', target_chain='A'):
        """
        sequence based alignment of target to reference
        :return:
        """
        target.detect_ligand_bonds()
        target.add_hydrogens()
        binding_site_superposition = Protein.ChainSuperposition()
        binding_site_superposition.superpose(reference[reference_chain],
                                             target[target_chain])

        return target


def main():
    results = Organiser()
    print np.mean(results.rmsd())

if __name__ == "__main__":
    main()