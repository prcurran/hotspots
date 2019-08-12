import os
import pickle
import shutil
import sys

import pandas as pd
from ccdc import io
from ccdc.cavity import Cavity
from ccdc.protein import Protein
from hotspots.atomic_hotspot_calculation import _AtomicHotspot, _AtomicHotspotResult
from hotspots.calculation import Runner, Buriedness, ExpBuriedness
from hotspots.grid_extension import Grid
from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.hs_utilities import Helper
from hotspots.result import Results, Extractor


def create_directory(path):
    """
    create a directory if it doesn't already exist

    :param path:
    :return: path
    """
    if not os.path.exists(path):
        os.mkdir(path)

    return path


class HotspotPipeline(object):

    def __init__(self, apo, buriedness_method, protein_id, ligand_id):
        """
        Initilising the HotspotPipeline object will set the structure of the output files. If an alternative
        naming scheme is required, change here.

        :param str apo: PDB code of the structure to be calculated
        :param str buriedness_method: either 'ghecom', 'ligsite', 'ghecom_internal' determines the buriesness method
        to be employed
        """
        # inputs
        self.apo = apo
        self.buriedness_method = buriedness_method
        self.protein_id = protein_id
        self.ligand_id = [l.split("_") for l in ligand_id]
        # outputs

        # directories
        self.working_dir_base_base = create_directory(os.path.join("data", self.apo[1:3]))
        self.working_dir_base = create_directory(os.path.join(self.working_dir_base_base, self.apo))
        self.working_dir = create_directory(os.path.join(self.working_dir_base, self.buriedness_method))

        # files

        # 'hotspot' protein files
        self.log_file = os.path.join(self.working_dir_base_base, "{}_{}.log".format(self.apo, self.buriedness_method))
        self.apo_protein = os.path.join(self.working_dir_base, "{}.pdb".format(self.apo))
        self.apo_prep = os.path.join(self.working_dir_base, "{}_prep.pdb".format(self.apo))
        self.buriedness = os.path.join(self.working_dir, "buriedness_{}.grd".format(self.buriedness_method))

        # 'other' protein files
        self.other_pdbs = {p: os.path.join(self.working_dir_base, "{}.pdb".format(p)) for p in self.protein_id}
        self.extracted_ligands = {pid: {j: os.path.join(self.working_dir_base, "{}.mol2".format(j))
                                        for j in self.ligand_id[i]}
                                  for i, pid in enumerate(self.protein_id)}

        self.ligand_volume = {pid: {j: os.path.join(self.working_dir_base_base, "{}.volume".format(j))
                                    for j in self.ligand_id[i]}
                              for i, pid in enumerate(self.protein_id)}

        self.cavity_rank = os.path.join(self.working_dir, "cavity.rank")

        # these get set during the run depending on the number of cavities detected
        self.cavities = []
        self.superstar = []
        self.hotspot = []
        self.bcv = {}
        self.all_overlaps = []

    def _prep_protein(self):
        """
        removes no structural ligands and solvents from the protein. Hotspot method requires the cavitiy to be empty

        :return: None
        """

        # input
        prot = Protein.from_file(self.apo_protein)

        # task
        prot.remove_all_waters()
        prot.detect_ligand_bonds()
        for l in prot.ligands:
            if 'HEM' not in l.identifier:
                prot.remove_ligand(l.identifier)

        # output
        with io.MoleculeWriter(self.apo_prep) as w:
            w.write(prot)

    def extract_ligands(self, other_id, lig_id):
        """
        extracts the relevant ligand(s) from the aligned PDB to a mol2 file

        :param str other_id: position in list of 'other' proteins
        :return:
        """
        # inputs
        other = Protein.from_file(self.other_pdbs[other_id])

        # tasks
        other.detect_ligand_bonds()
        print([a.identifier for a in other.ligands])
        print(other_id, lig_id)
        relevant = [l for l in other.ligands if lig_id in l.identifier.split(":")[1][0:3]]

        # output
        with io.MoleculeWriter(self.extracted_ligands[other_id][lig_id]) as writer:
            writer.write(relevant[0])                               # if more than one ligand detected, output the first

    def _get_ligand_volume(self, other_id, lig_id):
        """
        from a ligand, output a molecular volume in A^3

        :param i: position in list of 'other' proteins
        :return:
        """
        # inputs
        ligand = io.MoleculeReader(self.extracted_ligands[other_id][lig_id])[0]

        # tasks
        g = Grid.from_molecule(ligand)
        vol = g.count_grid() * (g.spacing ** 3)

        # output
        with open(self.ligand_volume[other_id][lig_id], 'w') as f:
            f.write(str(vol))

    def _get_buriedness_grid(self):
        """
        calculates the buriedness grid

        (if self.buriedness_method = ligsite, do nothing as we get this grid for free during the SuperStar Calc)
        :return: None
        """

        # inputs
        prot = Protein.from_file(self.apo_prep)

        # tasks
        coords = [a.coordinates for a in prot.atoms]
        out_grid = Grid.initalise_grid(coords=coords, padding=2, spacing=0.5)

        if self.buriedness_method == 'ghecom':
            b = Buriedness(protein=prot, out_grid=out_grid)
            g = b.calculate().grid
            shutil.rmtree(b.settings.working_directory)

        elif self.buriedness_method == 'ghecom_internal':
            b = ExpBuriedness(prot=prot, out_grid=out_grid)
            g = b.buriedness_grid()

        elif self.buriedness_method == 'ligsite':
            g = None
            pass

        else:
            raise TypeError("Not a valid pocket detection method")

        # outputs
        if self.buriedness_method == 'ligsite':
            pass
        else:
            if not os.path.exists(os.path.dirname(self.buriedness)):
                os.mkdir(os.path.dirname(self.buriedness))

            g.write(self.buriedness)

    def _get_cavities(self, min_vol):
        """
        detect cavities using Cavity API, generate new directory for each cavity

        :return: None
        """

        # inputs
        cavs = [c for c in Cavity.from_pdb_file(self.apo_prep) if c.volume > min_vol]

        # task, output
        cav_volume_dic = {os.path.join(self.working_dir, 'cavity_{}'.format(i), "cavity.volume"): c.volume
                   for i, c in enumerate(cavs)}

        cav_dic = {os.path.join(self.working_dir, 'cavity_{}'.format(i)): Helper.cavity_centroid(c)
                   for i, c in enumerate(cavs)}

        for path, vol in cav_volume_dic.items():
            with open(os.path.join(path), 'w') as f:
                f.write(str(vol))

        for path, origin in cav_dic.items():
            create_directory(path)
            with open(os.path.join(path, "cavity_origin.pkl"), 'wb') as handle:
                pickle.dump(origin, handle)

        # update attr
        self.cavities_volumes = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "cavity.volume") for p in
                                range(len(cav_dic))]
        self.cavities = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "cavity_origin.pkl") for p in range(len(cav_dic))]
        self.superstar = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "superstar") for p in range(len(cav_dic))]
        self.hotspot = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "hotspot") for p in range(len(cav_dic))]
        self.cavity_score = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "cavity.score") for p in range(len(cav_dic))]

        self.bcv = {i: {pid: {k: os.path.join(self.working_dir, "cavity_{}".format(i), "bcv",  "volume_{}".format(k))
                              for k in self.ligand_id[j]}
                        for j, pid in enumerate(self.protein_id)}
                    for i in range(len(cav_dic))}
        self.all_overlaps = {i: {pid: {k: os.path.join(self.working_dir, "cavity_{}".format(i), "bcv",  "overlap_{}.percentage".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in range(len(cav_dic))}

        self.atomic_overlaps = {i: {pid: {k: os.path.join(self.working_dir, "cavity_{}".format(i), "bcv",  "atomic_overlap_{}.dat".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in range(len(cav_dic))}

        self.matched = {i: {pid: {k: os.path.join(self.working_dir, "cavity_{}".format(i), "bcv",  "atom_match_{}.percentage".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in range(len(cav_dic))}

    def _get_superstar(self, cav_id=None):
        """
        calculate SuperStar for each cavity

        if the buriedness method is ligsite, write out the grid for later

        :param cav_id:
        :return:
        """
        # input
        prot = Protein.from_file(self.apo_prep)
        with open(self.cavities[cav_id], 'rb') as handle:
            cavity_origin = pickle.load(handle)

        # tasks
        a = _AtomicHotspot()
        a.settings.atomic_probes = {"apolar": "AROMATIC CH CARBON",
                                    "donor": "UNCHARGED NH NITROGEN",
                                    "acceptor": "CARBONYL OXYGEN"}

        self.superstar_grids = a.calculate(prot, nthreads=None, cavity_origins=[cavity_origin])

        sr = Results(protein=prot,
                     super_grids={result.identifier: result.grid for result in self.superstar_grids})

        #  outputs
        if not os.path.exists(self.superstar[cav_id]):
            os.mkdir(self.superstar[cav_id])

        with HotspotWriter(path=self.superstar[cav_id], zip_results=True) as w:
            w.write(sr)

        if self.buriedness_method == 'ligsite':
            # only write if it doesn't exist i.e. the first cavity run
            if not os.path.exists(self.buriedness):
                for ss in self.superstar_grids:
                    if ss.identifier == "apolar":
                        ss.buriedness.write(self.buriedness)

    def _get_hotspot(self, cav_id):
        """
        calculate hotspot map from pre-calculated superstar and buriedness grids

        :param cav_id:
        :return:
        """
        # inputs
        prot = Protein.from_file(self.apo_prep)
        sr = HotspotReader(path=os.path.join(self.superstar[cav_id], "out.zip")).read()
        superstar = [_AtomicHotspotResult(identifier=ident, grid=grid, buriedness=None)
                     for ident, grid in sr.super_grids.items()]
        buriedness = Grid.from_file(self.buriedness)

        # tasks
        h = Runner()

        s = h.Settings()
        s.apolar_translation_threshold = 14
        s.polar_translation_threshold = 14
        s.polar_contributions = False
        s.sphere_maps = False
        s.nrotations = 3000

        hr = h.from_superstar(prot, superstar, buriedness, settings=s, clear_tmp=True)

        # output
        if not os.path.exists(self.hotspot[cav_id]):
            os.mkdir(self.hotspot[cav_id])

        with HotspotWriter(self.hotspot[cav_id], zip_results=True) as writer:
            writer.write(hr)

    def _get_bcv(self, cav_id, other_id, lig_id):
        """
        generate a BCV for each cavity, and each required volume

        :param cav_id:
        :return:
        """
        # inputs
        hr = HotspotReader(path=os.path.join(self.hotspot[cav_id], "out.zip")).read()
        with open(self.ligand_volume[other_id][lig_id], 'r') as f:
            target_volume = f.read()

        # task
        extractor = Extractor(hr)
        bcv = extractor.extract_volume(volume=int(float(target_volume)))

        # output
        out = self.bcv[cav_id][other_id][lig_id]

        create_directory(os.path.dirname(out))
        create_directory(out)

        with HotspotWriter(path=out, grid_extension=".grd", zip_results=True) as writer:
            writer.write(bcv)

    def _get_volume_overlap(self, cav_id, other_id, lig_id):
        """
        find the highest median bcv from all cavities, calculate percentage over between the best bcv
        and each query ligand

        :return:
        """
        # inputs
        mol = io.MoleculeReader(self.extracted_ligands[other_id][lig_id])[0]
        path = os.path.join(self.bcv[cav_id][other_id][lig_id], "out.zip")
        if os.path.exists(path):
            hr = HotspotReader(path).read()

            # tasks
            sg = Grid.get_single_grid(hr.super_grids, mask=False)
            other = Grid.from_molecule(mol)
            overlap = sg.percentage_overlap(other=other)

            # output
            with open(self.all_overlaps[cav_id][other_id][lig_id], 'w') as writer:
                writer.write(str(overlap))
        else:
            print("no BCV for cavity {}, BCV {}".format(cav_id, lig_id))

    def _get_atomic_overlap(self, cav_id, other_id, lig_id):
        """
        find the highest median bcv from all cavities, calculate percentage over between the best bcv
        and each query ligand

        :return:
        """
        # inputs
        mol = io.MoleculeReader(self.extracted_ligands[other_id][lig_id])[0]
        path = os.path.join(self.bcv[cav_id][other_id][lig_id], "out.zip")
        if os.path.exists(path):
            hr = HotspotReader(path).read()

            # tasks
            out = hr.atomic_volume_overlap(mol)

        else:
            print("no BCV for cavity {}, BCV {}".format(cav_id, lig_id))
            out = {"donor": {}, "acceptor": {}, "apolar": {}}
            for a in mol.heavy_atoms:
                t = Helper.get_atom_type(a)
                if t == "doneptor":
                    out["donor"].update({a.label: 0.0})
                    out["acceptor"].update({a.label: 0.0})
                else:
                    out[t].update({a.label: 0.0})

        # output
        with open(self.atomic_overlaps[cav_id][other_id][lig_id], 'w') as writer:
            writer.write(str(out))

    def _get_matched_atoms(self, cav_id, other_id, lig_id):
        """


        :param cav_id:
        :param other_id:
        :param lig_id:
        :return:
        """
        # inputs
        mol = io.MoleculeReader(self.extracted_ligands[other_id][lig_id])[0]
        path = os.path.join(self.bcv[cav_id][other_id][lig_id], "out.zip")
        if os.path.exists(path):
            hr = HotspotReader(path).read()

            # tasks
            perc, type_dic = hr.percentage_matched_atoms(mol=mol, threshold=25, match_atom_types=True)

            # output
            with open(self.matched[cav_id][other_id][lig_id], 'w') as writer:
                writer.write(str(perc) + "\n")
                writer.write(str(type_dic))
        else:
            print("no BCV for cavity {}, BCV {}".format(cav_id, lig_id))

    def _score_cavity(self, cav_id):
        """
        score the cavity using the hotspot score

        :param cav_id:
        :return:
        """
        print(self.apo)
        hr = HotspotReader(os.path.join(self.hotspot[cav_id], "out.zip")).read()
        s = hr.score()

        with open(self.cavity_score[cav_id], "w") as f:
            f.write(str(s))

    def _rank_cavities(self):
        """
        rank the best continuous volumes by ligands

        :return:
        """
        # inputs
        cavity_by_score = {}
        for i in range(len(self.cavity_score)):
            with open(self.cavity_score[i], "r") as f:
                cavity_by_score.update({f.read():i})

        # task
        top = sorted(cavity_by_score.keys(), reverse=True)[0]

        # output
        with open(self.cavity_rank, "w") as writer:
            writer.write(str(cavity_by_score[top]))

    def run(self, rerun=False):
        # step 1: prepare protein
        self._prep_protein()

        # step 2: extract ligands from "other" proteins
        for other_id, lig_dic in self.extracted_ligands.items():
            for lig_id, path in lig_dic.items():
                if not os.path.exists(path) or rerun:
                    self.extract_ligands(other_id=other_id, lig_id=lig_id)

        # step 3: get extracted ligand volumes
                if not os.path.exists(self.ligand_volume[other_id][lig_id]) or rerun:
                    self._get_ligand_volume(other_id=other_id, lig_id=lig_id)

        # step 4: calculate pocket buriedness
        if not os.path.exists(self.buriedness) or rerun:
            self._get_buriedness_grid()

        # step 5: cavity detection
        self._get_cavities(min_vol=200)

        # step 6: superstar calculation
        for cav_id, prot_dic in self.bcv.items():
            if not os.path.exists(os.path.join(self.superstar[cav_id], "out.zip")) or rerun:
                self._get_superstar(cav_id=cav_id)

        # step 7: hotspot calculation
            if not os.path.exists(os.path.join(self.hotspot[cav_id], "out.zip")) or rerun:
                self._get_hotspot(cav_id=cav_id)

        # step 7a: score cavities
            if not os.path.exists(self.cavity_score[cav_id]):
                self._score_cavity(cav_id=cav_id)

        # step 8: bcv calculatuion
            for prot_id, lig_dic in prot_dic.items():
                for lig_id, path in lig_dic.items():
                    if not os.path.exists(self.bcv[cav_id][prot_id][lig_id]) or rerun:
                        try:
                            self._get_bcv(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)
                        except:
                            continue

        # step 9: overlap analysis
                    if not os.path.exists(self.all_overlaps[cav_id][prot_id][lig_id]) or rerun:
                        self._get_volume_overlap(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)

        # step 9a: atomic overlap analysis
                    print("Atomic Overlap {}".format(self.apo))
                    if not os.path.exists(self.atomic_overlaps[cav_id][prot_id][lig_id] or rerun):
                        self._get_atomic_overlap(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)

        # step 10: atom match analysis
                    if not os.path.exists(self.matched[cav_id][prot_id][lig_id]) or rerun:
                        self._get_matched_atoms(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)

        # step 11: get cavity rank
        if not os.path.exists(self.cavity_rank):
            self._rank_cavities()


def main():
    # inputs
    buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']

    for method in buriedness_methods:

        # try:
        #
        #     hp = HotspotPipeline(
        #                          apo=sys.argv[1],
        #                          buriedness_method=method,
        #                          protein_id=sys.argv[2].split(","),
        #                          ligand_id=sys.argv[3].split(",")
        #     )
        #
        #     hp.run(rerun=False)
        #
        # except:
        #     continue


        hp = HotspotPipeline(
                             apo=sys.argv[1],
                             buriedness_method=method,
                             protein_id=sys.argv[2].split(","),
                             ligand_id=sys.argv[3].split(",")
        )

        hp.run(rerun=False)


if __name__ == '__main__':
    main()
