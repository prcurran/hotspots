from __future__ import division

import os
import pickle
import shutil
import sys
import time

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
        self.apo_protein = os.path.join(self.working_dir_base, "{}.pdb".format(self.apo))
        self.apo_prep = os.path.join(self.working_dir_base, "{}_prep.pdb".format(self.apo))
        self.buriedness = os.path.join(self.working_dir, "buriedness_{}.grd".format(self.buriedness_method))
        self.buriedness_time = os.path.join(self.working_dir, "buriedness_{}.time".format(self.buriedness_method))

        # 'other' protein files
        self.other_pdbs = {p: os.path.join(self.working_dir_base, "{}.pdb".format(p)) for p in self.protein_id}
        self.extracted_ligands = {pid: {j: os.path.join(self.working_dir_base, "{}.mol2".format(j))
                                        for j in self.ligand_id[i]}
                                  for i, pid in enumerate(self.protein_id)}
        self.ligand_volume = {pid: {j: os.path.join(self.working_dir_base_base, "{}.volume".format(j))
                                    for j in self.ligand_id[i]}
                              for i, pid in enumerate(self.protein_id)}

        # analysis files
        self.cavity_rank = os.path.join(self.working_dir, "cavity.rank")
        self.ligand_cavity = os.path.join(self.working_dir, "ligand_cavity.dat")

        # these get set during the run depending on the number of cavities detected
        self.runs = ["global"]     # + cavities
        create_directory(os.path.join(self.working_dir, "global"))

        self.cavities = []
        self.cavities_volumes = []
        self.cavities = {}
        self.bounding_box = []

        self.superstar = {}
        self.hotspot = {}
        self.bcv = {}
        self.superstar_time = {}
        self.hotspot_time = {}
        self.bcv_time = {}
        self.bcv_threshold = {}

        self.bcv_lig_overlaps = {}
        self.bcv_hot_overlaps = {}
        self.hot_lig_overlaps = {}
        self.hot_hot_overlaps = {}

        self.matched = {}

    # CALCULATION

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
        start = time.time()
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
        finish = time.time()

        # outputs
        with open(self.buriedness_time, 'w') as t:
            t.write(str(finish - start))

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

        # task
        for i in range(len(cavs)):
            create_directory(path=os.path.join(self.working_dir, 'cavity_{}'.format(i)))

        cav_dic = {os.path.join(self.working_dir, 'cavity_{}'.format(i)): Helper.cavity_centroid(c)
                   for i, c in enumerate(cavs)}

        cav_volume_dic = {os.path.join(self.working_dir, 'cavity_{}'.format(i), "cavity.volume"): c.volume
                   for i, c in enumerate(cavs)}

        cav_bb = {os.path.join(self.working_dir, 'cavity_{}'.format(i), "bounding_box.pkl"): c.bounding_box
                   for i, c in enumerate(cavs)}

        # output
        for path, origin in cav_dic.items():
            with open(os.path.join(path, "cavity_origin.pkl"), 'wb') as handle:
                pickle.dump(origin, handle)

        for path, vol in cav_volume_dic.items():
            with open(os.path.join(path), 'w') as f:
                f.write(str(vol))

        for path, bb in cav_bb.items():
            with open(os.path.join(path), 'wb') as h:
                pickle.dump(bb, h)

        # update attr
        self.runs += ["cavity_{}".format(i) for i in range(len(cavs))]

        self.cavities = {"cavity_{}".format(p): os.path.join(self.working_dir, 'cavity_{}'.format(p), "cavity_origin.pkl")
                         for p in range(len(cav_dic))}
        self.cavities_volumes = {"cavity_{}".format(p): os.path.join(self.working_dir, 'cavity_{}'.format(p), "cavity.volume")
                                 for p in range(len(cav_dic))}
        self.cavity_score = {"cavity_{}".format(p): os.path.join(self.working_dir, "cavity_{}".format(p), "cavity.score")
                             for p in range(len(cav_dic))}
        self.bounding_box = {"cavity_{}".format(p): os.path.join(self.working_dir, 'cavity_{}'.format(p), "bounding_box.pkl")
                             for p in range(len(cav_dic))}

        self.superstar = {p: os.path.join(self.working_dir, p, "superstar") for p in self.runs}
        self.superstar_time = {k: os.path.join(v, "time.time") for k, v in self.superstar.items()}

        self.hotspot = {p: os.path.join(self.working_dir, p, "hotspot") for p in self.runs}
        self.hotspot_time = {k: os.path.join(v, "time.time") for k, v in self.hotspot.items()}

        self.bcv = {i: {pid: {k: os.path.join(self.working_dir, i, "bcv",  "volume_{}".format(k))
                              for k in self.ligand_id[j]}
                        for j, pid in enumerate(self.protein_id)}
                    for i in self.runs}
        self.bcv_time = {i: {pid: {k: os.path.join(self.working_dir, i, "bcv",  "volume_{}", "time.time".format(k))
                                   for k in self.ligand_id[j]}
                             for j, pid in enumerate(self.protein_id)}
                         for i in self.runs}
        self.bcv_threshold = {i: {pid: {k: os.path.join(self.working_dir, i, "bcv",  "volume_{}, threshold.dat".format(k))
                                        for k in self.ligand_id[j]}
                                  for j, pid in enumerate(self.protein_id)}
                             for i in self.runs}

        self.bcv_lig_overlaps = {i: {pid: {k: os.path.join(self.working_dir, i, "bcv",  "lig_overlap_{}.percentage".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in self.runs}

        self.bcv_hot_overlaps = {i: {pid: {k: os.path.join(self.working_dir, i, "bcv",  "hot_overlap_{}.percentage".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in self.runs}

        self.hot_lig_overlaps = {i: {pid: {k: os.path.join(self.working_dir, i, "hotspot",  "lig_overlap_{}.percentage".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in self.runs}

        self.hot_hot_overlaps = {i: {pid: {k: os.path.join(self.working_dir, i, "hotspot",  "hot_overlap_{}.percentage".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in self.runs}

        self.matched = {i: {pid: {k: os.path.join(self.working_dir, i, "bcv",  "atom_match_{}.percentage".format(k))
                                       for k in self.ligand_id[j]}
                                 for j, pid in enumerate(self.protein_id)}
                             for i in self.runs}

    def _get_superstar(self, cav_id=None):
        """
        calculate SuperStar for each cavity

        if the buriedness method is ligsite, write out the grid for later

        :param cav_id:
        :return:
        """
        # input

        prot = Protein.from_file(self.apo_prep)

        if cav_id is 'global':
            cavity_origin = None
        else:
            with open(self.cavities[cav_id], 'rb') as handle:
                cavity_origin = [pickle.load(handle)]

        # tasks
        start = time.time()
        a = _AtomicHotspot()
        a.settings.atomic_probes = {"apolar": "AROMATIC CH CARBON",
                                    "donor": "UNCHARGED NH NITROGEN",
                                    "acceptor": "CARBONYL OXYGEN"}

        self.superstar_grids = a.calculate(prot, nthreads=None, cavity_origins=cavity_origin)

        sr = Results(protein=prot,
                     super_grids={result.identifier: result.grid for result in self.superstar_grids})
        finish = time.time()

        #  outputs
        if not os.path.exists(self.superstar[cav_id]):
            os.mkdir(self.superstar[cav_id])

        if cavity_id is not 'global':
            out = os.path.join(a.settings.temp_dir, str(0))
        else:
            out = a.settings.temp_dir

        for interaction in ["apolar", "acceptor", "donor"]:
            shutil.copyfile(os.path.join(out, "{}.cavity.mol2".format(interaction)),
                            os.path.join(self.superstar[cav_id], "{}.cavity.mol2".format(interaction)))

        shutil.make_archive(os.path.join(self.superstar[cav_id], "superstar"), 'zip', out)

        with HotspotWriter(path=self.superstar[cav_id], zip_results=True) as w:
            w.write(sr)

        with open(self.superstar_time[cav_id], 'w') as t:
            t.write(str(finish - start))

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
        start = time.time()
        h = Runner()

        s = h.Settings()
        s.apolar_translation_threshold = 14
        s.polar_translation_threshold = 14
        s.polar_contributions = False
        s.sphere_maps = False
        s.nrotations = 3000

        hr = h.from_superstar(prot, superstar, buriedness, settings=s, clear_tmp=True)
        finish = time.time()
        # output
        if not os.path.exists(self.hotspot[cav_id]):
            os.mkdir(self.hotspot[cav_id])

        with open(self.hotspot_time[cav_id], 'w') as t:
            t.write(str(finish - start))

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
        start = time.time()
        extractor = Extractor(hr)
        bcv = extractor.extract_volume(volume=int(float(target_volume)))
        finish = time.time()

        # output
        out = self.bcv[cav_id][other_id][lig_id]

        create_directory(os.path.dirname(out))
        create_directory(out)

        with HotspotWriter(path=out, grid_extension=".grd", zip_results=True) as writer:
            writer.write(bcv)

        with open(self.bcv_time[cav_id][other_id][lig_id], 'w') as t:
            t.write(str(finish - start))

        with open(self.bcv_threshold[cav_id][other_id][lig_id], 'w') as s:
            s.write(str(bcv.step_threshold))

    # ANALYSIS

    def _get_volume_overlap(self, cav_id, other_id, lig_id):
        """
        find the highest median bcv from all cavities, calculate percentage over between the best bcv
        and each query ligand

        :return:
        """
        def nonzero(val):
            if val == 0:
                return 1
            else:
                return val
        # inputs
        mol = io.MoleculeReader(self.extracted_ligands[other_id][lig_id])[0]
        path1 = os.path.join(self.hotspot[cav_id], "out.zip")
        path2 = os.path.join(self.bcv[cav_id][other_id][lig_id], "out.zip")
        thresholds = [10, 14, 17]

        if os.path.exists(path1) and os.path.exists(path2):
            bcv = HotspotReader(path2).read()
            hot = HotspotReader(path1).read()

            # tasks
            other = Grid.from_molecule(mol)

            bcv_sg = Grid.get_single_grid(bcv.super_grids, mask=False)
            bcv_overlap = bcv_sg._mutually_inclusive(other=other).count_grid()

            lig_vol = (other > 0).count_grid()
            bcv_vol = (bcv_sg > 0).count_grid()

            hot_sgs = [(Grid.get_single_grid(hot.super_grids, mask=False) > t)
                       for t in thresholds]
            hot_vols = [nonzero(hot_sg.count_grid())
                        for hot_sg in hot_sgs]
            hot_overlap = [hot_sg._mutually_inclusive(other=other).count_grid() for hot_sg in hot_sgs]

            # output
            with open(self.bcv_lig_overlaps[cav_id][other_id][lig_id], 'w') as writer:
                writer.write(str((bcv_overlap / lig_vol) * 100))

            with open(self.bcv_hot_overlaps[cav_id][other_id][lig_id], 'w') as writer:
                writer.write(str((bcv_overlap / bcv_vol) * 100))

            with open(self.hot_lig_overlaps[cav_id][other_id][lig_id], 'w') as writer:
                hot_lig = [str((a / lig_vol) * 100) for a in hot_overlap]
                print hot_lig
                writer.write(",".join(hot_lig))

            with open(self.hot_hot_overlaps[cav_id][other_id][lig_id], 'w') as writer:
                hot_hot = [str((hot_overlap[i] / hot_vols[i]) * 100) for i in range(len(thresholds))]
                writer.write(",".join(hot_hot))

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
        This is the ligand overlap implimentation in the DoGsiter paper

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
            perc, type_dic = hr.percentage_matched_atoms(mol=mol, threshold=0, match_atom_types=True)

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

    def _get_ligand_cavity(self):
        """
        determine which cavity contains the ligand, output cav_id to file

        :return:
        """
        lc = []
        tolerance = 0
        mols = [MoleculeReader(path)[0] for other_id, lig_dic in self.extracted_ligands.items()
                for l, path in lig_dic.items()]

        point = [round(np.mean([a.coordinates.x for mol in mols for a in mol.heavy_atoms])),
                 round(np.mean([a.coordinates.y for mol in mols for a in mol.heavy_atoms])),
                 round(np.mean([a.coordinates.z for mol in mols for a in mol.heavy_atoms]))]

        for i in range(len(self.bounding_box)):
            with open(self.bounding_box[i], 'rb') as handle:
                bb = pickle.load(handle)
            
            mini = bb[0]
            maxi = bb[1]
            
            if all([mini.x - tolerance < point[0] < maxi.x + tolerance,
                    mini.y - tolerance < point[1] < maxi.y + tolerance,
                    mini.z - tolerance < point[2] < maxi.z + tolerance]):

                lc.append(str(i))

        with open(self.ligand_cavity, "w") as f:
            f.write(",".join(lc))

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
        for cav_id, prot_dic in self.bcv.items():       # cav_id is run_id
            if not os.path.exists(os.path.join(self.superstar[cav_id], "out.zip")) or rerun:
                self._get_superstar(cav_id=cav_id)

        # step 7: hotspot calculation
            if not os.path.exists(os.path.join(self.hotspot[cav_id], "out.zip")) or rerun:
                self._get_hotspot(cav_id=cav_id)

        # step 7a: score cavities
            if cav_id != 'global':
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
                    if not os.path.exists(self.hot_hot_overlaps[cav_id][prot_id][lig_id]) or rerun:
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
    # buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']

    buriedness_methods = ['ghecom']

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

        hp.run(rerun=True)


if __name__ == '__main__':
    main()
