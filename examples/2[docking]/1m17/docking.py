"""
This is an example script to accompany our communication on Hotspot API.
    - Use Case 1: Hotspot-Guided Docking

Author: P. R. Curran
"""
from __future__ import print_function
import argparse
import os
import tempfile

import numpy as np

from ccdc.cavity import Cavity
from ccdc.docking import Docker
from ccdc.io import MoleculeWriter, _CSDDatabaseLocator
from ccdc.molecule import Molecule, Atom
from ccdc.protein import Protein
from ccdc.utilities import _private_importer
with _private_importer():
    import ChemicalAnalysisLib
    import ConformerGeneratorLib

from hotspots import calculation
from hotspots import hs_io
from hotspots import hs_docking
from hotspots import hs_utilities
from hotspots import result
from hotspots.pdb_python_api import _Ligand, PDBResult


class Organiser(argparse.ArgumentParser):
    """
    class to handle docking run with hotspot insights
    """

    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        # handle command line arguments
        self.add_argument(
            'path',
            help='path to working directory'
        )

        self.add_argument(
            'pdb',
            help='PDB code for target'
        )

        self.add_argument(
            'chemical_id',
            help='PDB code for target'
        )

        self.add_argument(
            '-hs', '--hotspot_guided',
            default=True,
            help='Use Hotspot insights to guide docking'
        )

        self.args = self.parse_args()

        # create temp for output files
        self.args.path = hs_io.Helper.get_out_dir(self.args.path)
        self.temp = tempfile.mkdtemp()

        # calculate hotspot using Hotspots API
        if self.args.hotspot_guided is True:
            try:
                self.hr = hs_io.HotspotReader(path=os.path.join(self.args.path, "out.zip")).read()

            except IOError:
                h = calculation.Runner()
                settings = h.Settings()
                settings.nrotations = 3000
                settings.sphere_maps = False
                self.hr = h.from_pdb(pdb_code=self.args.pdb,
                                     charged_probes=False,
                                     probe_size=3,
                                     buriedness_method='ghecom',
                                     nprocesses=5,
                                     settings=settings)

                with hs_io.HotspotWriter(path=os.path.join(self.args.path), zip_results=True) as hw:
                    hw.write(self.hr)

        # generate molecule for docking
        self.search_ligands = os.path.join(self.temp, self.args.chemical_id + ".mol2")
        self.ligand = self.from_smiles(smiles=_Ligand.from_chemicalid(chemicalid=self.args.chemical_id).smiles,
                                       path=self.search_ligands,
                                       identifier=self.args.chemical_id)

        # dock search ligands into hotspot protein
        self.docked_ligands = self.dock()

        if self.args.hotspot_guided is True:
            self.rescored_ligands = self.rescore()

    @staticmethod
    def from_smiles(smiles, path, identifier=None, generate_initial_sites=True):
        """
        Create a :class:`ccdc.molecule.Molecule` from a SMILES string.

        *e.g.*::

             ethene = Molecule.from_smiles('C=C', 'Ethene')

        If ``identifier`` is not specified, the SMILES string will be used as the
        molecule identifier.

        :param smiles: str
        :param identifier: str
        :param generate_initial_sites: boolean - whether to include an initial guess at 3D coordinates
        :return: a :class:`ccdc.molecule.Molecule` instance with coordinates
        """
        if identifier is None:
            identifier = smiles

        if generate_initial_sites:
            parameter_files = _CSDDatabaseLocator.get_conformer_parameter_file_location()
            molmaker = ConformerGeneratorLib.MoleculeTo3D(parameter_files)
            mol = Molecule(identifier, molmaker.create_conformation(smiles))
        else:
            molmaker = ChemicalAnalysisLib.SMILESMoleculeMaker()
            mol = Molecule(identifier, _molecule=molmaker.siteless_atoms(smiles))

        with MoleculeWriter(path) as w:
            w.write(mol)

        return mol

    def rescore(self):
        """
        Rescores docking poses in Fragment Hotspot Map fields.
        :return: list of :class:`ccdc.molecule.Molecule`
        """
        rescored_mols = []
        mol_dic = {}
        self.docked_ligands = [self.hr.score(m.molecule, tolerance=0) for m in self.docked_ligands]
        for m in self.docked_ligands:
            array = []
            for atm in m.atoms:
                if atm.atomic_number > 0:
                    if type(atm.partial_charge) is not float:
                        array.append(0)
                    else:
                        array.append(atm.partial_charge)
            if len(array) == 0:
                continue

            else:
                g = np.mean(np.array(array))
                try:
                    mol_dic[g].append(m)
                except KeyError:
                    mol_dic.update({g: [m]})

        for key in sorted(mol_dic.keys(), reverse=True):
            for mol in mol_dic[key]:
                mol.identifier = "{}".format(key)
                rescored_mols.append(mol)


        return rescored_mols

    def dock(self):
        """
        Setup and execution of docking run with GOLD.

        NB: Docking Settings class is imported from the Hotspots API rather than Docking API. This is essential for
        running hotspot guided docking.
        :return: a :class:`ccdc.io.MoleculeReader`
        """
        docker = Docker()
        docker.settings = hs_docking.DockerSettings()

        # download protein
        PDBResult(self.args.pdb).download(self.temp)
        protein = Protein.from_file(os.path.join(self.temp, self.args.pdb + ".pdb"))
        protein.remove_all_waters()
        protein.remove_all_metals()
        protein.add_hydrogens()
        for l in protein.ligands:
            protein.remove_ligand(l.identifier)

        f = os.path.join(self.temp, self.args.pdb + ".mol2")
        with MoleculeWriter(f) as w:
            w.write(protein)

        # setup
        docker.settings.add_protein_file(f)

        # create binding site from list of residues
        cavs = Cavity.from_pdb_file(os.path.join(self.temp, self.args.pdb + ".pdb"))
        cavs[0].to_pymol_file("test.py")
        c = {}
        for i, cav in enumerate(cavs):
            cav.feats = []
            for f in cav.features:
                try:
                    cav.feats.append(f.residue)
                except:
                    continue

            # cav.feats = [f.residue for f in cav.features]
            cav.len = len(cav.feats)
            c.update({cav.len: cav.feats})
            cav.to_pymol_file("{}.py".format(i))


        selected_cavity = max(c.keys())

        docker.settings.binding_site = docker.settings.BindingSiteFromListOfResidues(protein=docker.settings.proteins[0],
                                                                                     residues=c[selected_cavity])
        docker.settings.fitness_function = 'plp'
        docker.settings.autoscale = 100.
        docker.settings.output_directory = self.temp
        docker.settings.output_file = "docked_ligands.mol2"
        docker.settings.add_ligand_file(self.search_ligands, ndocks=25)

        # constraints
        if self.args.hotspot_guided is True:
            e_settings = result.Extractor.Settings()
            e_settings.mvon = True
            extractor = result.Extractor(self.hr, settings=e_settings)
            bv = extractor.extract_best_volume(volume=300)[0]
            f = hs_utilities.Helper.get_out_dir(os.path.join(self.args.path, "best_volume"))

            with hs_io.HotspotWriter(path=f) as hw:
                hw.write(bv)

            constraints = docker.settings.HotspotHBondConstraint.create(protein=docker.settings.proteins[0],
                                                                        hr=bv,
                                                                        weight=5,
                                                                        min_hbond_score=0.2,
                                                                        max_constraints=5)

            for constraint in constraints:
                docker.settings.add_constraint(constraint)
            docker.settings.add_fitting_points(hr=bv)

            mol = Molecule(identifier="constraints")
            for constraint in constraints:
                for a in constraint.atoms:
                    mol.add_atom(Atom(atomic_symbol="C",
                                      atomic_number=14,
                                      label="Du",
                                      coordinates=a.coordinates))

            with MoleculeWriter(os.path.join(self.args.path, "constaints.mol2")) as w:
                w.write(mol)

        docker.dock()
        results = docker.Results(docker.settings)
        return results.ligands

    def write(self, fname="results.mol2"):
        """
        Writes out docking pose.
        :param fname: str, path to output file
        :return: None
        """
        with MoleculeWriter(os.path.join(self.args.path, fname)) as w:
            try:
                for l in self.rescored_ligands:
                    w.write(l)
            except:
                for l in self.docked_ligands:
                    mol = l.molecule
                    mol.identifier = "{}".format(l.fitness())
                    w.write(mol)

    def clean_up(self):
        """
        clean up files generated from docking
        :return:
        """
        os.remove(self.args.pdb + ".mol2")

        try:
            os.rename("fit_pts.mol2", os.path.join(self.args.path, "fit_pts.mol2"))
        except OSError:
            pass

        os.rename("api_gold.conf", os.path.join(self.args.path, "api_gold.conf"))

        try:
            os.rename("cavity.residues", os.path.join(self.args.path, "cavity.residues"))
        except OSError:
            pass

def main():
    result = Organiser()
    result.write()
    result.clean_up()


if __name__ == "__main__":
    main()
