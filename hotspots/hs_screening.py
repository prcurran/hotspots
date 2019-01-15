"""
The :mod:`hs_screening` module contains classes which
faciliate the applications of insights from Fragment
Hotspot Maps to ccdc's virtual screening tools

The main classes of the :mod:`hs_screening` module are:

- :class:``
- :class:``

"""
import tempfile
import os

from ccdc import docking
from ccdc import protein
from ccdc import molecule
from ccdc import io
from ccdc.utilities import _private_importer

with _private_importer():
    import DockingLib

DockingLib.licence_check()

from best_volume import Extractor


class DockerSettings(docking.Docker.Settings):
    class HotspotHBondConstraint(docking.Docker.Settings.Constraint):
        """A protein HBond constraint constructed from a hotspot
        Assign Protein Hbond constraints based on the highest scoring interactions."""

        def __init__(self, atoms, weight=5.0, min_hbond_score=0.001, hotspot_score=None, _constraint=None):
            '''Initialise a HotspotHBondConstraint constraint.

            :param atoms: a list :class:`ccdc.molecule.Atom` instances from the protein.
             The atoms should be donatable hydrogens or acceptor atoms.
            :param weight: the penalty to be applied for no atom of the list forming an HBond.
            :param min_hbond_score: the minimal score of an HBond to be considered a valid HBond.
            '''
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
            self.hotspot_scores = hotspot_score
            self.weight = weight
            self.min_hbond_score = min_hbond_score
            if _constraint is None:
                _constraint = DockingLib.GoldProteinHBondConstraint()
                _constraint.from_string(self._to_string())
            super(self.__class__, self).__init__(_constraint)

        @staticmethod
        def from_hotspot(protein, hr, max_constraints=2, weight=5.0, min_hbond_score=0.001):
            """
            construct a Hotspot HBondConstraint
            :param hr: a `hotspots.calculation.Result`
            :return:
            """
            for atm in protein.atoms:
                atm.partial_charge = 0

            prot = hr.score(protein)
            with io.MoleculeWriter("/home/pcurran/use_case/scored_prot.mol2") as w:
                w.write(prot)
            atm_dic = {atm.partial_charge: atm for atm in prot.atoms
                       if type(atm.partial_charge) is float
                       and ((atm.atomic_number == 1 and atm.neighbours[0].is_donor) or atm.is_acceptor)
                       }

            scores = sorted([f[0] for f in atm_dic.items()], reverse=True)[:max_constraints]

            selected = [atm_dic[s] for s in scores]
            return DockerSettings.HotspotHBondConstraint(atoms=selected,
                                                         weight=weight,
                                                         hotspot_score=scores,
                                                         min_hbond_score=min_hbond_score)

        def _to_string(self):
            s = '%.4f %.4f %s' % (
                self.weight, self.min_hbond_score, ' '.join(str(a.index + 1) for a in self.atoms)
            )
            return s

    def add_apolar_fitting_points(self, hr, volume=400, threshold=17, mode='threshold'):
        """
        help on volumes
        sets apolar fitting points for docking run
        :param hr: `hotspots.calculation.Result`
        :param volume: int
        :param threshold: int
        :param mode: str, 'threshold' or 'bcv'
        :return:
        """
        temp = tempfile.mkdtemp()
        mol = molecule.Molecule(identifier="fitting_pts")

        if mode == 'threshold':
            dic = hr.super_grids["apolar"].grid_value_by_coordinates(threshold=17)

        elif mode == 'bcv':
            extractor = Extractor(hr=hr, volume=volume, mode="global")
            bv = extractor.extracted_hotspots[0].best_island
            dic = bv.grid_value_by_coordinates(threshold=17)

        else:
            raise TypeError("{} not supported, see documentation for detials".format(mode))

        for score, v in dic.items():
            for pts in v:
                atm = molecule.Atom(atomic_symbol='C',
                                    atomic_number=14,
                                    label='{:.2f}'.format(score),
                                    coordinates=pts)
                atm.partial_charge = score
                mol.add_atom(atom=atm)

        fname = os.path.join(temp, "test.mol2")
        with io.MoleculeWriter(fname) as w:
            w.write(mol)

        self.fitting_points_file = fname


docking.Docker.Settings = DockerSettings
