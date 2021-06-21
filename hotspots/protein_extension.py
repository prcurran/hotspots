from ccdc import protein
import numpy as np
from scipy.spatial import distance


def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length


class Protein(protein.Protein):
    def __init__(self, identifier, _molecule=None, _protein_structure=None):
        super().__init__(identifier, _molecule=_molecule, _protein_structure=_protein_structure)

    class BindingSiteFromGrid(protein.Protein.BindingSiteFromListOfAtoms):
        def __init__(self, prot, grid, within=4, background_val=0):
            list_of_residues = self._detect_from_grid(prot, grid, within, background_val)
            self.list_of_atoms = tuple(a for r in list_of_residues for a in r.atoms)
            self.list_of_residues = list_of_residues
            super().__init__(prot, self.list_of_atoms)

        @staticmethod
        def _detect_from_grid(prot, grid, within, background_val=0):
            """
            detects a binding site from a grid object

            :param prot: the protein from which the binding site is to be extracted
            :param grid: a grid, for best results this should be a mask

            :type prot: `ccdc.protein.Protein`
            :type grid: `hotspots.grid_extension.Grid`

            :return: list of residues
            :rtype: list
            """
            # threshold set to 3 as this is the approx cavity def cutoff
            residue_centroids = [centroid(np.array([atm.coordinates for atm in res.atoms])) for res in prot.residues]
            residues = prot.residues
            bs_residues = set()

            for pnt in grid.edge_detection(edge_definition=background_val):
                dists = distance.cdist(np.array([pnt]), np.array(residue_centroids))[0]
                print(dists)
                index_array = np.nonzero(dists < within)[0]
                if not len(index_array) == 0:
                    for i in index_array:
                        bs_residues.add(residues[i])
            return list(bs_residues)


protein.Protein = Protein


class Residue(protein.Protein.Residue):
    def __init__(self, i, _residue):
        super().__init__(i, _residue)

    def is_incomplete(self):
        # alternate atom location may be described - use set
        aa_atm_dic = {'MET': 8, 'GLU': 9, 'ASN': 8, 'PHE': 11, 'GLN': 9, 'LYS': 9,
                      'VAL': 7, 'ILE': 8, 'GLY': 4, 'THR': 7, 'TYR': 12, 'ALA': 5,
                      'ARG': 11, 'LEU': 8, 'SER': 6, 'HIS': 10, 'PRO': 7, 'ASP': 8,
                      'CYS': 6, 'TRP': 14}
        res_id = self.identifier.split(":")[1][:3]
        if res_id not in aa_atm_dic:
            return False
        else:
            if len({atm.label for atm in self.atoms}) < aa_atm_dic[res_id]:
                return True
            else:
                return False


protein.Protein.Residue = Residue