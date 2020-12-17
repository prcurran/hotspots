import numpy as np
from scipy.spatial import distance

from ccdc import protein

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
        def __init__(self, prot, grid, within=4):
            list_of_residues = self._detect_from_grid(prot, grid, within)
            self.list_of_atoms = tuple(a for r in list_of_residues for a in r.atoms)
            self.list_of_residues = list_of_residues
            super().__init__(prot, self.list_of_atoms)

        @staticmethod
        def _detect_from_grid(prot, grid, within):
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

            for pnt in grid.edge_detection():
                dists = distance.cdist(np.array([pnt]), np.array(residue_centroids))[0]
                index_array = np.nonzero(dists < within)[0]
                if not len(index_array) == 0:
                    for i in index_array:
                        bs_residues.add(residues[i])
            return list(bs_residues)


protein.Protein = Protein