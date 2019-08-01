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
        top = sorted(cavity_by_score.keys(), reverse=True)

        # obj_dic = {}
        # for cav_id, prot_dic in self.bcv.items():
        #     for prot_id, lig_dic in prot_dic.items():
        #         for lig_id, path in lig_dic.items():
        #             if os.path.exists(os.path.join(path, "out.zip")):
        #                 if cav_id not in obj_dic:
        #                     obj_dic.update({cav_id: {}})
        #                 obj_dic[cav_id].update({lig_id: HotspotReader(os.path.join(path, "out.zip")).read()})
        #
        # # tasks
        # lines = []
        # cavs = obj_dic.keys()
        # ligands = set([lig for cav in cavs for lig in obj_dic[cav].keys()])
        # print(obj_dic)
        #
        # for ligand in ligands:
        #     cav_by_score = {obj_dic[cav_id][ligand].score(): cav_id for cav_id in cavs if ligand in obj_dic[cav_id]}
        #     top_cavity = cav_by_score[sorted(cav_by_score.keys(), reverse=True)[0]]
        #     lines.append("{}: {}\n".format(ligand, top_cavity))

        # output
        with open(self.cavity_rank, "w") as writer:
            writer.write(str(top))

    def _log_message(self, message=""):
        m = "{}\n".format(message)
        with open(self.log_file, 'a') as f:
            f.write(m)

    def run(self, rerun=False):
        p = ",".join(self.protein_id)
        l = ",".join([l for prot in self.ligand_id for l in prot])
        lines = ["#\n",
                 "PDB file used: {}\n".format(self.apo),
                 "Buriedness method used: {}\n".format(self.buriedness_method),
                 "ligands {} from {} respectively\n".format(l, p),
                 "#\n"
                 ]

        with open(self.log_file, 'w+') as f:
            f.writelines(lines)

        # step 1: prepare protein
        self._prep_protein()
        self._log_message("Remove ligands: passed")

        # step 2: extract ligands from "other" proteins
        for other_id, lig_dic in self.extracted_ligands.items():
            for lig_id, path in lig_dic.items():
                if not os.path.exists(path) or rerun:
                    self.extract_ligands(other_id=other_id, lig_id=lig_id)
                    self._log_message("Ligand Extracted: passed, ID = {}".format(lig_id))

        # step 3: get extracted ligand volumes
        for other_id, lig_dic in self.ligand_volume.items():
            for lig_id, path in lig_dic.items():
                if not os.path.exists(path) or rerun:
                    self._get_ligand_volume(other_id=other_id, lig_id=lig_id)
                    self._log_message("Ligand Volume: passed, ID = {}".format(lig_id))

        # step 4: calculate pocket buriedness
        if not os.path.exists(self.buriedness) or rerun:
            self._get_buriedness_grid()
            self._log_message("Buriedness calculation: passed")

        # step 5: cavity detection
        self._get_cavities(min_vol=200)
        self._log_message("Cavity calculation: passed")

        # step 6: superstar calculation
        for i, path in enumerate(self.superstar):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                self._get_superstar(cav_id=i)
                self._log_message("SuperStar calculation: passed, ID = cavity_{}".format(i))

        # step 7: hotspot calculation
        for i, path in enumerate(self.hotspot):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                self._get_hotspot(cav_id=i)
                self._log_message("Hotspot calculation: passed, ID = cavity_{}".format(i))

        # step 8: bcv calculatuion
        for cav_id, prot_dic in self.bcv.items():
            for prot_id, lig_dic in prot_dic.items():
                for lig_id, path in lig_dic.items():
                    if not os.path.exists(path) or rerun:
                        try:
                            self._get_bcv(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)
                            self._log_message("BCV calculation: passed, CAV_ID = {}, LIG_ID = {}".format(cav_id, lig_id))
                        except:
                            self._log_message(
                                "BCV calculation: failed, CAV_ID = {}, LIG_ID = {}".format(cav_id, lig_id))

        # step 9: overlap analysis
        for cav_id, prot_dic in self.all_overlaps.items():
            for prot_id, lig_dic in prot_dic.items():
                for lig_id, path in lig_dic.items():
                    if not os.path.exists(path) or rerun:
                        self._get_volume_overlap(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)
                        self._log_message("Volume overlap analysis: passed")

        # step 10: atom match analysis
        for cav_id, prot_dic in self.matched.items():
            for prot_id, lig_dic in prot_dic.items():
                for lig_id, path in lig_dic.items():
                    if not os.path.exists(path) or rerun:
                        self._get_matched_atoms(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)
                        self._log_message("Match atom analysis: passed")

        # step 11: rank cavities
        for i, path in enumerate(self.cavity_score):
            if not os.path.exists(path):
                self._score_cavity(cav_id=i)

        if not os.path.exists(self.cavity_rank):
            self._rank_cavities()

    def _cav_volume(self, cav):
        if os.path.exists(self.cavities_volumes[cav]):
            with open(self.cavities_volumes[cav], "r") as f:
                vol = f.read()
            return float(vol)
        else:
            return "NaN"

    def _step_1_status(self):
        if os.path.exists(self.apo_prep):
            return 1
        else:
            return 0

    def _step_2_status(self, other_id, lig_id):
        if os.path.exists(self.extracted_ligands[other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_3_status(self, other_id, lig_id):
        if os.path.exists(self.ligand_volume[other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_4_status(self):
        if os.path.exists(self.buriedness):
            return 1
        else:
            return 0

    def _step_5_status(self, cav):
        if os.path.exists(self.cavities[cav]):
            return 1
        else:
            return 0

    def _step_6_status(self, cav):
        if os.path.exists(os.path.join(self.superstar[cav], "out.zip")):
            return 1
        else:
            return 0

    def _step_7_status(self, cav):
        if os.path.exists(os.path.join(self.hotspot[cav], "out.zip")):
            return 1
        else:
            return 0

    def _step_8_status(self, cav, other_id, lig_id):
        if os.path.exists(os.path.join(self.bcv[cav][other_id][lig_id], "out.zip")):
            return 1
        else:
            return 0

    def _step_9_status(self, cav, other_id, lig_id):
        if os.path.exists(self.all_overlaps[cav][other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_10_status(self, cav, other_id, lig_id):
        if os.path.exists(self.matched[cav][other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_11_status(self):
        if os.path.exists(self.cavity_rank):
            return 1
        else:
            return 0

    def status(self):
        """
        report the status of the calculations

        :return:
        """
        keys = ["apo", "buriedness_method", "cavity_id", "other_id", "ligand_id", "cav_vol", "apo_prep", "extract_ligand",
                "extract_vol", "calculate_buriedness", "cavity_detection", "calculate_superstar", "calculate_hotspots", "calculate_bcv",
                "volume_overlap", "matched_atoms", "cavity_rank"]
        data = []
        self._get_cavities(min_vol=200)

        s1 = self._step_1_status()
        s4 = self._step_4_status()
        s11 = self._step_11_status()

        for cav in range(len(self.cavities)):
            vol = self._cav_volume(cav)
            s5 = self._step_5_status(cav)
            s6 = self._step_6_status(cav)
            s7 = self._step_7_status(cav)
            for i, prot_id in enumerate(self.protein_id):
                for lig_id in self.ligand_id[i]:
                    s2 = self._step_2_status(prot_id, lig_id)
                    s3 = self._step_3_status(prot_id, lig_id)
                    s8 = self._step_8_status(cav, prot_id, lig_id)
                    s9 = self._step_9_status(cav, prot_id, lig_id)
                    s10 = self._step_10_status(cav, prot_id, lig_id)

                    data.append([self.apo, self.buriedness_method, cav, prot_id, lig_id, vol, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11])

        return pd.DataFrame(dict(zip(keys, (zip(*data)))))

    def _get_vol_data(self, cav, other_id, lig_id):
        """
        extract the volume data

        :param int cav: cavity identifier
        :param str other_id: protein identifier
        :param str lig_id: ligand identifier
        :return:
        """
        if os.path.exists(self.all_overlaps[cav][other_id][lig_id]):
            with open(self.all_overlaps[cav][other_id][lig_id], "r") as f:
                vol = f.read()
            return float(vol)
        else:
            return 0.0

    def _is_top_cavity(self, cav):
        """
        determine if top ranked cavity

        :param cav:
        :return:
        """
        if os.path.exists(self.cavity_rank):
            with open(self.cavity_rank, "r") as f:
                c = dict([y for y in [x.split(":") for x in f.read().split("\n")] if len(y) > 1])
                concensus = list(set(c.values()))
                if len(concensus) > 1:
                    print(self.apo)
                    print(concensus)
                    raise IOError('Take a look Pedro')

                else:
                    if int(cav) == int(concensus[0]):
                        return True
                    else:
                        return False
            return 1
        else:
            return "NAN"

    def analysis(self):
        """
        report the volume analysis

        :return:
        """
        keys = ["apo", "buriedness_method", "cavity_id", "other_id", "ligand_id", "bcv_calc", "top_cavity", "volume_overlap"]
        data = []
        self._get_cavities(min_vol=200)

        for cav in range(len(self.cavities)):
            v1 = self._is_top_cavity(cav)
            for i, prot_id in enumerate(self.protein_id):
                for lig_id in self.ligand_id[i]:
                    s8 = self._step_8_status(cav, prot_id, lig_id)

                    v2 = self._get_vol_data(cav, prot_id, lig_id)
                    data.append([self.apo, self.buriedness_method, cav, prot_id, lig_id, s8, v1, v2])

        return pd.DataFrame(dict(zip(keys, (zip(*data)))))


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
