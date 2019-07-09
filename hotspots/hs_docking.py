"""
The :mod:`hotspots.hs_docking` module contains functionality which
faciliates the **automatic** application of insights from Fragment
Hotspot Maps to docking.


This module is designed to extend the existing CSD python API

More information about the CSD python API is available:
    - The Cambridge Structural Database C.R. Groom, I. J. Bruno, M. P. Lightfoot and S. C. Ward, Acta Crystallographica Section B, B72, 171-179, 2016 [DOI: 10.1107/S2052520616003954]
    - CSD python API 2.0.0 `documentation <https://downloads.ccdc.cam.ac.uk/documentation/API/>`_

More information about the GOLD method is available:
    - Development and Validation of a Genetic Algorithm for Flexible Docking G. Jones, P. Willett, R. C. Glen, A. R. Leach and R. Taylor, J. Mol. Biol., 267, 727-748, 1997 [DOI: 10.1006/jmbi.1996.0897]
"""
from __future__ import print_function
import os
import tempfile

import numpy as np
from ccdc import docking
from ccdc import io
from ccdc import molecule
from ccdc.utilities import _private_importer
from hotspots.result import Extractor, Results

with _private_importer():
    import DockingLib

DockingLib.licence_check()


class DockerSettings(docking.Docker.Settings):
    """
    A class to handle the integration of Fragment Hotspot Map data with GOLD

    This class is designed to mirror the existing CSD python API for smooth integration. For use, import this class
    as the docking settings rather than directly from the Docking API.

    >>> from ccdc.docking import Docker
    >>> from ccdc.protein import Protein
    >>> from hotspots.calculation import Runner
    >>> from hotspots.hs_docking import DockerSettings

    >>> protein = Protein.from_file("1hcl.pdb")

    >>> runner = Runner()
    >>> hs = runner.from_protein(protein)

    >>> docker.settings.add_protein_file("1hcl.pdb")
    >>> docker.settings.add_ligand_file("dock_me.mol2", ndocks=25)
    >>> constraints = docker.settings.HotspotHBondConstraint.from_hotspot(protein=docker.settings.proteins[0], hr=hs)
    >>> docker.settings.add_constraint(constraints)
    >>> docker.dock()

    >>> docker.Results(docker.settings).ligands

    """
    class HotspotHBondConstraint(docking.Docker.Settings.Constraint):
        """
        A protein HBond constraint constructed from a hotspot
        Assign Protein Hbond constraints based on the highest scoring interactions.

        :param list atoms: list of :class:`ccdc.molecule.Atom` instances from the protein. *NB: The atoms should be donatable hydrogens or acceptor atoms.*
        :param weight: the penalty to be applied for no atom of the list forming an HBond.
        :param min_hbond_score: the minimal score of an HBond to be considered a valid HBond.
        """
        def __init__(self, atoms, weight=5.0, min_hbond_score=0.001, _constraint=None):
            self.atoms = atoms[:]
            for a in self.atoms:
                if not self._is_protein_atom(a):
                    raise RuntimeError('One of the atoms is not in the protein')
                if (
                        not (a.atomic_number == 1 and a.neighbours[0].is_donor) and
                        not (a.is_acceptor)
                ):
                    raise RuntimeError('One of the atoms does not form an HBond')
            self._add_to_protein = atoms[0]._molecule
            self.weight = weight
            self.min_hbond_score = min_hbond_score
            if _constraint is None:
                _constraint = DockingLib.GoldProteinHBondConstraint()
                _constraint.from_string(self._to_string())
            super(self.__class__, self).__init__(_constraint)

        @staticmethod
        def _from_file(path, protein, weight, min_hbond_score=0.001, max=2):
            """
            create a hotspot constraint from file

            :return:
            """
            constraints = Results._ConstraintData.read(path)

            return [DockerSettings.HotspotHBondConstraint(atoms=[protein.atoms[index]],
                                                          weight=float(weight),
                                                          min_hbond_score=float(min_hbond_score))
                    for i, index in enumerate(constraints.score_by_index.values()) if i < max]

        @staticmethod
        def create(protein, hr, max_constraints=2, weight=5.0, min_hbond_score=0.001, cutoff=10):
            """
            creates a :class:`hotspots.hs_docking.HotspotHBondConstraint`

            :param `ccdc.protein.Protein` protein: the protein to be used for docking
            :param `hotspots.calculation.Result` hr: a result from Fragment Hotspot Maps
            :param int max_constraints: max number of constraints
            :param float weight: the constraint weight (default to be determined)
            :param float min_hbond_score: float between 0.0 (bad) and 1.0 (good) determining the minimum hydrogen bond quality in the solutions.
            :param cutoff: minimum score required to assign the constraint
            :return list: list of :class:`hotspots.hs_docking.HotspotHBondConstraint`
            """
            # for atm in protein.atoms:
            #     atm.partial_charge = int(0)
            # prot = hr.score(protein)
            #
            # coords = np.array([a.coordinates for a in prot.atoms])
            # atm_dic = {atm.partial_charge: atm for atm in prot.atoms
            #            if type(atm.partial_charge) is float
            #            and ((atm.atomic_number == 1 and atm.neighbours[0].is_donor) or atm.is_acceptor)
            #            and _is_solvent_accessible(coords, atm, min_distance=2)
            #            }
            #
            # print(atm_dic)
            #
            # if len(atm_dic) > max_constraints:
            #     scores = sorted([f[0] for f in atm_dic.items()], reverse=True)[:max_constraints]
            # else:
            #     scores = sorted([f[0] for f in atm_dic.items()], reverse=True)

            constraints = hr._docking_constraint_atoms(p=protein, max_constraints=max_constraints)

            return [DockerSettings.HotspotHBondConstraint(atoms=[protein.atoms[constraints.score_by_index[a]]],
                                                          weight=weight,
                                                          min_hbond_score=min_hbond_score)
                    for a in constraints.score_by_index.keys() if a > cutoff]

        def _to_string(self):
            s = '%.4f %.4f %s' % (
                self.weight, self.min_hbond_score, ' '.join(str(a.index + 1) for a in self.atoms)
            )
            return s

    def _add_fitting_points_from_file(self, path='fit_pts.mol2'):
        """
        creates fitting pts from a .mol2 file

        :param str path: path to fitting pts file
        """
        self.fitting_points_file = path

    def generate_fitting_points(self, hr, volume=400, threshold=17, mode='threshold'):
        """
        uses the Fragment Hotspot Maps to generate GOLD fitting points.

        GOLD fitting points are used to help place the molecules into the protein cavity. Pre-generating these fitting
        points using the Fragment Hotspot Maps helps to biast results towards making Hotspot interactions.

        :param `hotspots.result.Result` hr: a Fragment Hotspot Maps result
        :param int volume: volume of the occupied by fitting points in Angstroms ^ 3
        :param float threshold: points above this value will be included in the fitting points
        :param str mode: 'threshold'- assigns fitting points based on a score cutoff or 'bcv'- assigns fitting points from best continuous volume analysis (recommended)
        """
        temp = tempfile.mkdtemp()
        mol = molecule.Molecule(identifier="fitting_pts")

        if mode == 'threshold':
            dic = hr.super_grids["apolar"].grid_value_by_coordinates(threshold=threshold)

        elif mode == 'bcv':
            extractor = Extractor(hr=hr, volume=volume, mode="global")
            bv = extractor.extracted_hotspots[0].best_island
            dic = bv.grid_value_by_coordinates(threshold=17)

        else:
            raise TypeError("{} not supported, see documentation for details".format(mode))

        for score, v in dic.items():
            for pts in v:
                atm = molecule.Atom(atomic_symbol='C',
                                    atomic_number=14,
                                    label='{:.2f}'.format(score),
                                    coordinates=pts)
                atm.partial_charge = score
                mol.add_atom(atom=atm)

        fname = os.path.join(temp, 'fit_pts.mol2')
        with io.MoleculeWriter(fname) as w:
            w.write(mol)

        self.fitting_points_file = fname


def _is_solvent_accessible(protein_coords, atm, min_distance=2):
    """
    given a protein and an atom of a protein, determine if the atom is accessible to solvent

    :param protein:
    :param atm:
    :return:
    """
    if str(atm.atomic_symbol) == 'H':
        atm_position = np.array(atm.coordinates)
        neighbour = np.array(atm.neighbours[0].coordinates)
        direction = np.subtract(atm_position, neighbour) * 2
        position = np.array([direction + atm_position])
        distance = min(np.linalg.norm(protein_coords - position, axis=1))
        if distance > min_distance:
            return True
        else:
            return False

    else:
        return True


docking.Docker.Settings = DockerSettings
