import os
import pickle
import shutil
import subprocess
import sys


PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
else:
    import urllib2

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

    def __init__(self, pdb, buriedness_method, protein_id, ligand_id):
        """
        Initilising the HotspotPipeline object will set the structure of the output files. If an alternative
        naming scheme is required, change here.

        :param str pdb: PDB code of the structure to be calculated
        :param str buriedness_method: either 'ghecom', 'ligsite', 'ghecom_internal' determines the buriesness method
        to be employed
        """
        # inputs
        self.pdb = pdb
        self.buriedness_method = buriedness_method
        self.protein_id = protein_id
        self.ligand_id = [l.split("_") for l in ligand_id]
        # outputs

        # directories
        self.working_dir_base_base = create_directory(os.path.join("pdb_files", self.pdb[1:3]))
        self.working_dir_base = create_directory(os.path.join("pdb_files", self.pdb[1:3], self.pdb))
        self.working_dir = create_directory(os.path.join("pdb_files", self.pdb[1:3], self.pdb, self.buriedness_method))

        # files

        # 'hotspot' protein files
        self.log_file = os.path.join(self.working_dir_base_base, "{}_{}.log".format(self.pdb, self.buriedness_method))
        self.biological_assembly = os.path.join(self.working_dir_base, "biological_assembly.pdb")
        self.protonated = os.path.join(self.working_dir_base, "protonated.pdb")
        self.no_wat = os.path.join(self.working_dir_base, "protonated_no_wat.pdb")
        self.buriedness = os.path.join(self.working_dir, "buriedness_{}.grd".format(self.buriedness_method))

        # 'other' protein files
        self.other_pdbs = {p: os.path.join(self.working_dir_base, "{}.pdb".format(p)) for p in self.protein_id}
        self.aligned_pdbs = {p: os.path.join(self.working_dir_base, "{}_aligned_to_{}.pdb".format(p, self.pdb))
                             for p in self.protein_id}

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

    def _download_pdb(self):
        """
        download the pdb from the RCSB. pdb1 indicates the protein is the biologial assembly

        :return: None
        """

        # task
        url = "https://files.rcsb.org/download/{}.pdb1".format(self.pdb)
        print(url)
        pdbfile = urllib2.urlopen(url).read()

        # output
        with open(self.biological_assembly, 'wb') as out_file:
            out_file.write(pdbfile)

    def _protonate(self):
        """
        protonate the downloaded biological assembley. (using pdb2pqr)

        :return: None
        """

        # input, task, output
        cmd = "/vagrant/pdb2pqr-linux-bin64-2.1.1/pdb2pqr --ff=amber --chain {} {}".format(self.biological_assembly,
                                                                                           self.protonated)

        subprocess.call(cmd, shell=sys.platform != 'win32')

        if not os.path.exists(self.protonated):
            raise RuntimeError("PDB2PQR protonation failed")

    def _protonate_backup(self):
        """
        if pdb2pqr fails, just use the CSD python API.
        (Sometimes, missing atoms on the end of loops cause pdb2pqr to fall over.)

        :return:
        """

        # input
        prot = Protein.from_file(self.biological_assembly)

        # task
        prot.add_hydrogens()

        # output
        with io.MoleculeWriter(self.protonated) as writer:
            writer.write(prot)

    def _remove_wat_lig(self):
        """
        removes no structural ligands and solvents from the protein. Hotspot method requires the cavitiy to be empty

        :return: None
        """

        # input
        prot = Protein.from_file(self.protonated)

        # task
        prot.remove_all_waters()
        prot.detect_ligand_bonds()
        for l in prot.ligands:
            if 'HEM' not in l.identifier:
                prot.remove_ligand(l.identifier)

        # output
        with io.MoleculeWriter(self.no_wat) as w:
            w.write(prot)

    def _download_other_pdb(self, other_id):
        """
        download the pdb from the RCSB. pdb1 indicates the protein is the biologial assembly

        :return: None
        """

        # task

        url = "https://files.rcsb.org/download/{}.pdb1".format(self.protein_id[other_id])
        print(url)
        pdbfile = urllib2.urlopen(url).read()

        # output
        with open(self.other_pdbs[other_id], 'wb') as out_file:
            out_file.write(pdbfile)

    def _align_other_to_main(self, other_id):
        """
        aligns 'other' proteins to 'hotspot' protein. alignment done per chain

        :param str other_id: position in list of 'other' proteins
        :return:
        """
        # input
        hotspot = Protein.from_file(self.no_wat)
        other = Protein.from_file(self.other_pdbs[other_id])

        # task
        hotspot_chain = [c.identifier for c in hotspot.chains][0]           # take first chain

        other.add_hydrogens()
        other.detect_ligand_bonds()

        relevant = [l for l in other.ligands
                    for a in self.ligand_id
                    for b in a
                    if b in l.identifier.split(":")[1][0:3]]                # align other PDB by relevant ligand's chain

        other_chain = relevant[0].identifier.split(":")[0]

        binding_site_superposition = Protein.ChainSuperposition()
        bs_rmsd, bs_transformation = binding_site_superposition.superpose(hotspot[hotspot_chain],
                                                                          other[other_chain])

        self._log_message("Alignment RMSD: {}".format(bs_rmsd))

        # output
        with io.MoleculeWriter(self.aligned_pdbs[other_id]) as writer:
            writer.write(other)

    def extract_ligands(self, other_id, lig_id):
        """
        extracts the relevant ligand(s) from the aligned PDB to a mol2 file

        :param str other_id: position in list of 'other' proteins
        :return:
        """
        # inputs
        other = Protein.from_file(self.aligned_pdbs[other_id])

        # tasks
        relevant = [l for l in other.ligands if lig_id in l.identifier.split(":")[1][0:3]]

        # output
        with io.MoleculeWriter(self.extracted_ligands[lig_id]) as writer:
            writer.write(relevant[0])                               # if more than one ligand detected, output the first

    def _get_ligand_volume(self, other_id, lig_id):
        """
        from a ligand, output a molecular volume in A^3

        :param i: position in list of 'other' proteins
        :return:
        """
        # inputs
        ligand = io.MoleculeReader(self.extract_ligands(other_id=other_id, lig_id=lig_id))[0]

        # tasks
        g = Grid.from_molecule(ligand)
        vol = g.count_grid() * (g.spacing ** 3)
        print(vol)

        # output
        with open(self.ligand_volume[lig_id], 'w') as f:
            f.write(str(vol))

    def _get_buriedness_grid(self):
        """
        calculates the buriedness grid

        (if self.buriedness_method = ligsite, do nothing as we get this grid for free during the SuperStar Calc)
        :return: None
        """

        # inputs
        prot = Protein.from_file(self.no_wat)

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
        cavs = [c for c in Cavity.from_pdb_file(self.no_wat) if c.volume > min_vol]

        # task, output
        cav_dic = {os.path.join(self.working_dir, 'cavity_{}'.format(i)): Helper.cavity_centroid(c)
                   for i, c in enumerate(cavs)}

        for path, origin in cav_dic.items():
            create_directory(path)

            with open(os.path.join(path, "cavity_origin.pkl"), 'wb') as handle:
                pickle.dump(origin, handle)

        # update attr
        self.cavities = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "cavity_origin.pkl") for p in range(len(cav_dic))]
        self.superstar = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "superstar") for p in range(len(cav_dic))]
        self.hotspot = [os.path.join(self.working_dir, 'cavity_{}'.format(p), "hotspot") for p in range(len(cav_dic))]
        self.bcv = {i: {k: os.path.join(self.working_dir, "cavity_{}".format(i), "bcv",  "volume_{}".format(k))
                        for j in self.ligand_id for k in j} for i in range(len(cav_dic))}
        self.all_overlaps = {i: {pid: {k: os.path.join(self.working_dir, "cavity_{}".format(i), "bcv",  "overlap_{}.percentage".format(k))
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
        prot = Protein.from_file(self.no_wat)
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
        prot = Protein.from_file(self.no_wat)
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

    def _get_bcv(self, cav_id, lig_id):
        """
        generate a BCV for each cavity, and each required volume

        :param cav_id:
        :return:
        """
        # inputs
        hr = HotspotReader(path=os.path.join(self.hotspot[cav_id], "out.zip")).read()
        with open(self.ligand_volume[lig_id], 'r') as f:
            target_volume = f.read()

        # task
        extractor = Extractor(hr)
        bcv = extractor.extract_volume(volume=int(float(target_volume)))

        # output
        out = self.bcv[cav_id][lig_id]

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
        path = os.path.join(self.bcv[cav_id][lig_id], "out.zip")
        if os.path.exists(path):
            hr = HotspotReader(path).read()

            # tasks
            sg = Grid.get_single_grid(hr.super_grids, mask=False)
            other = Grid.from_molecule(mol)
            overlap = sg.percentage_overlap(other=other)
            print(overlap)

            # output
            with open(self.all_overlaps[cav_id][other_id][lig_id], 'w') as writer:
                writer.write(str(overlap))
        else:
            print("no BCV for cavity {}, BCV {}".format(cav_id, lig_id))

    def _rank_cavities(self):
        """
        rank the best continuous volumes by ligands

        :return:
        """
        # inputs
        obj_dic = {}
        for cav_id, lig_dic in self.bcv.items():
            for lig_id, path in lig_dic.items():
                if os.path.exists(os.path.join(path, "out.zip")):
                    if cav_id not in obj_dic:
                        obj_dic.update({cav_id: {}})
                    obj_dic[cav_id].update({lig_id: HotspotReader(os.path.join(path, "out.zip")).read()})

        print(obj_dic)

        # tasks
        lines = []
        cavs = obj_dic.keys()
        ligands = set([lig for cav in cavs for lig in obj_dic[cav].keys()])

        for ligand in ligands:
            cav_by_score = {obj_dic[cav_id][ligand].score(): cav_id for cav_id in cavs}
            top_cavity = cav_by_score[sorted(cav_by_score.keys(), reverse=True)[0]]
            lines.append("{}: {}\n".format(lig_id, top_cavity))

        # output
        with open(self.cavity_rank, "w") as writer:
            writer.writelines(lines)

    def _log_message(self, message=""):
        m = "{}\n".format(message)
        with open(self.log_file, 'a') as f:
            f.write(m)

    def run(self, rerun=False):
        p = ",".join(self.protein_id)
        l = ",".join([l for prot in self.ligand_id for l in prot])
        lines = ["#\n",
                 "PDB file used: {}\n".format(self.pdb),
                 "Buriedness method used: {}\n".format(self.buriedness_method),
                 "ligands {} from {} respectively\n".format(l, p),
                 "#\n"
                 ]

        with open(self.log_file, 'w+') as f:
            f.writelines(lines)

        try:
            if not os.path.exists(self.biological_assembly) or rerun:
                self._download_pdb()
            self._log_message("PDB Download: passed, ID = {}".format(self.pdb))
        except:
            self._log_message("PDB Download: failed, ID = {}".format(self.pdb))

        try:
            if not os.path.exists(self.protonated) or rerun:
                self._protonate()
            self._log_message("Protonation_1: passed")
        except:
            self._log_message("Protonation_1: failed")

            try:
                self._protonate_backup()
                self._log_message("Protonation_2: passed")
            except:
                self._log_message("Protonation_2: failed")

        try:
            if not os.path.exists(self.no_wat) or rerun:
                self._remove_wat_lig()
            self._log_message("Remove ligands: passed")

        except:
            self._log_message("Remove ligands: passed")

        for other_id, path in self.other_pdbs.items():
            if not os.path.exists(path) or rerun:
                try:
                    self._download_other_pdb(other_id)
                    self._log_message("PDB download: passed, ID = {}".format(other_id))
                except:
                    self._log_message("PDB download: failed, ID = {}".format(other_id))

        for other_id, path in self.aligned_pdbs.items():
            if not os.path.exists(path) or rerun:
                try:
                    self._align_other_to_main(other_id)
                    self._log_message("Alignment: passed, ID = {}".format(other_id))
                except:
                    self._log_message("PDB download: failed, ID = {}".format(other_id))

        for other_id, lig_dic in self.extracted_ligands.items():
            for lig_id, path in lig_dic.items():
                if not os.path.exists(path) or rerun:
                    try:
                        self.extract_ligands(other_id=other_id, lig_id=lig_id)
                        self._log_message("Ligand Extracted: passed, ID = {}".format(lig_id))
                    except:
                        self._log_message("Ligand Extracted: failed, ID = {}".format(lig_id))

        for other_id, lig_dic in self.ligand_volume.items():
            for lig_id, path in lig_dic.items():
                if not os.path.exists(path) or rerun:
                    try:
                        self._get_ligand_volume(other_id=other_id, lig_id=lig_id)
                        self._log_message("Ligand Volume: passed, ID = {}".format(lig_id))
                    except:
                        self._log_message("Ligand Volume: failed, ID = {}".format(lig_id))

        if not os.path.exists(self.buriedness) or rerun:
            try:
                self._get_buriedness_grid()
                self._log_message("Buriedness calculation: passed")
            except:
                self._log_message("Buriedness calculation: failed")

        try:
            self._get_cavities(min_vol=200)
            self._log_message("Cavity calculation: passed")
        except:
            self._log_message("Cavity calculation: failed")

        for i, path in enumerate(self.superstar):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                try:
                    self._get_superstar(cav_id=i)
                    self._log_message("SuperStar calculation: passed, ID = cavity_{}".format(i))
                except:
                    self._log_message("SuperStar calculation: failed, ID = cavity_{}".format(i))

        for i, path in enumerate(self.hotspot):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                try:
                    self._get_hotspot(cav_id=i)
                    self._log_message("Hotspot calculation: passed, ID = cavity_{}".format(i))
                except:
                    self._log_message("SuperStar calculation: failed, ID = cavity_{}".format(i))

        for cav_id, lig_dic in self.bcv.items():
            for lig_id, path in lig_dic.items():
                if not os.path.exists(path) or rerun:
                    try:
                        self._get_bcv(cav_id=cav_id, lig_id=lig_id)
                        self._log_message("BCV calculation: passed, CAV_ID = {}, LIG_ID = {}".format(cav_id, lig_id))
                    except:
                        self._log_message("BCV calculation: failed, CAV_ID = {}, LIG_ID = {}".format(cav_id, lig_id))

        for cav_id, prot_dic in self.all_overlaps.items():
            for prot_id, lig_dic in prot_dic.items():
                for lig_id, path in lig_dic.items():
                    if not os.path.exists(path) or rerun:
                        try:
                            self._get_volume_overlap(cav_id=cav_id, other_id=prot_id, lig_id=lig_id)
                            self._log_message("Volume overlap calculation: passed")
                        except:
                            self._log_message("Volume overlap calculation: failed")

        if not os.path.exists(self.cavity_rank):
            self._rank_cavities()

def main():
    # inputs
    df = pd.read_csv("readme.csv", index_col=0)
    buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']

    for index, row in df.iterrows():
        print(row.apo, index)
        for bm in buriedness_methods:
            hp = HotspotPipeline(pdb=row.apo,
                                 buriedness_method=bm,
                                 protein_id=[row.fragment, row.lead],
                                 ligand_id=[row.fragment_ID, row.lead_ID]
            )

            hp.run(rerun=False)
            if index >= 2:
                break


def test():
    # hp = HotspotPipeline(pdb='4J8N',
    #                      buriedness_method='ligsite',
    #                      protein_id=['2W1D', '2W1G'],
    #                      ligand_id=['L0D', 'L0G'])

    buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']
    for bm in buriedness_methods:
        hp = HotspotPipeline(pdb='3bqd',
                             buriedness_method=bm,
                             protein_id=['3BQD'],
                             ligand_id=['DAY'])

        hp.run(rerun=False)


if __name__ == '__main__':
    # main()
    test()