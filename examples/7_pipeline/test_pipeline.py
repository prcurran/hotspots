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
        self.working_dir_base_base = create_directory("pdb_files/{}".format(self.pdb[1:3]))
        self.working_dir_base = create_directory("pdb_files/{}/{}".format(self.pdb[1:3], self.pdb))
        self.working_dir = create_directory("pdb_files/{}/{}/{}".format(self.pdb[1:3], self.pdb,
                                                                        self.buriedness_method))

        # files

        # 'hotspot' protein files
        self.biological_assembly = "{}/biological_assembly.pdb".format(self.working_dir_base)
        self.protonated = "{}/protonated.pdb".format(self.working_dir_base)
        self.no_wat = "{}/protonated_no_wat.pdb".format(self.working_dir_base)
        self.buriedness = "{}/buriedness_{}.grd".format(self.working_dir, self.buriedness_method)

        # 'other' protein files
        self.other_pdbs = ["{}/{}.pdb".format(self.working_dir_base, p) for p in self.protein_id]
        self.aligned_pdbs = ["{}/{}_aligned_to_{}.pdb".format(self.working_dir_base, p, self.pdb)
                             for p in self.protein_id]

        self.extracted_ligands = [["{}/{}_from_{}.mol2".format(self.working_dir_base, self.protein_id[i], j)
                                   for j in self.ligand_id[i]]
                                  for i in range(len(self.protein_id))]

        self.volumes = [["{}/{}_from_{}.volume".format(self.working_dir_base, self.protein_id[i], j)
                         for j in self.ligand_id[i]]
                        for i in range(len(self.protein_id))]

        self.overlap = ["{}/{}_volume_overlap.percentage".format(self.working_dir, lig_id)
                        for p in self.ligand_id for lig_id in p]

        # these get set during the run depending on the number of cavities detected
        self.cavities = []
        self.superstar = []
        self.hotspot = []
        self.bcv = []
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

    def _download_other_pdb(self, i):
        """
        download the pdb from the RCSB. pdb1 indicates the protein is the biologial assembly

        :return: None
        """

        # task

        url = "https://files.rcsb.org/download/{}.pdb1".format(self.protein_id[i])
        print(url)
        pdbfile = urllib2.urlopen(url).read()

        # output
        with open(self.other_pdbs[i], 'wb') as out_file:
            out_file.write(pdbfile)

    def _align_other_to_main(self, i):
        """
        aligns 'other' proteins to 'hotspot' protein. alignment done per chain

        :param i: position in list of 'other' proteins
        :return:
        """
        # input
        hotspot = Protein.from_file(self.no_wat)
        other = Protein.from_file(self.other_pdbs[i])
        ligs = self.ligand_id[i]
        print(ligs)

        # task
        hotspot_chain = [c.identifier for c in hotspot.chains][0]       # this might cause problems

        other.add_hydrogens()
        other.detect_ligand_bonds()
        print(other.ligands)
        for l in other.ligands:
            print(l.identifier)
        relevant = [l for l in other.ligands if ligs[0] in l.identifier]
        print(relevant)                                         # multiple binders of interest will be in same chain
        other_chain = relevant[0].identifier.split(":")[0]

        binding_site_superposition = Protein.ChainSuperposition()
        bs_rmsd, bs_transformation = binding_site_superposition.superpose(hotspot[hotspot_chain],
                                                                          other[other_chain])

        print('RMSD: ', bs_rmsd)

        # output
        with io.MoleculeWriter(self.aligned_pdbs[i]) as writer:
            writer.write(other)

        mod_ligs = ["{}:{}".format(other_chain, lig) for lig in ligs]
        out_ligs = [l for l in other.ligands for mod_lig in mod_ligs if mod_lig in l.identifier]
        if len(out_ligs) > 1:
            print("More than 1 extracted ligand")

        for mol in out_ligs:
            for path in self.extracted_ligands[i]:
                if mol.identifier.split(":")[1][0:3] in path:
                    with io.MoleculeWriter(path) as writer:
                        writer.write(mol)

    def _get_ligand_volume(self, i):
        """
        from a ligand, output a molecular volume in A^3

        :param i: position in list of 'other' proteins
        :return:
        """

        # inputs
        ligands = [io.MoleculeReader(path)[0] for path in self.extracted_ligands[i]]

        # tasks
        for ligand in ligands:
            g = Grid.from_molecule(ligand)
            vol = g.count_grid() * (g.spacing ** 3)

            # output
            for path in self.volumes[i]:
                if ligand.identifier.split(":")[1][0:3] in path:
                    with open(path, 'wb') as f:
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
        self.cavities = [os.path.join(p, "cavity_origin.pkl") for p in cav_dic.keys()]
        self.superstar = [os.path.join(p, "superstar") for p in cav_dic.keys()]
        self.hotspot = [os.path.join(p, "hotspot") for p in cav_dic.keys()]

        vols = []
        for a, prot in enumerate(self.volumes):
            for b, ligand in enumerate(prot):
                with open(ligand, 'rb') as f:
                    vols.append((float(f.read()), self.ligand_id[a][b]))

        self.bcv_dir = [os.path.join(p, "bcv") for p in cav_dic.keys()]

        self.bcv = [[os.path.join(p, "bcv", "volume_{}".format(info[1])) for info in vols]
                    for p in cav_dic.keys()]

        self.all_overlaps = [os.path.join(b, "volume_overlap.percentage") for cav in self.bcv for b in cav]

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

        self.superstar_grids = a.calculate(prot,
                                           nthreads=None,
                                           cavity_origins=[cavity_origin])

        sr = Results(protein=prot,
                     super_grids={result.identifier: result.grid for result in self.superstar_grids}
                     )

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

    def _get_bcv(self, cav_id):
        """
        generate a BCV for each cavity, and each required volume

        :param cav_id:
        :return:
        """
        # inputs
        hr = HotspotReader(path=os.path.join(self.hotspot[cav_id], "out.zip")).read()
        vols = {}
        for a, prot in enumerate(self.volumes):
            for b, ligand in enumerate(prot):
                with open(ligand, 'rb') as f:
                    vols.update({float(f.read()): self.ligand_id[a][b]})

        # tasks
        for vol, ident in vols.items():
            try:
                e = Extractor(hr)
                bcv = e.extract_volume(volume=int(float(vol)))

                # output
                for path in self.bcv[cav_id]:
                    if ident in path:

                        for b in self.bcv_dir:
                            if not os.path.exists(b):
                                os.mkdir(b)

                        if not os.path.exists(path):
                            os.mkdir(path)

                        with HotspotWriter(path=path, grid_extension=".grd", zip_results=True) as writer:
                            writer.write(bcv)
            except:
                continue

    def _get_volume_overlap(self):
        """
        find the highest median bcv from all cavities, calculate percentage over between the best bcv
        and each query ligand

        :return:
        """

        # inputs
        mols = [io.MoleculeReader(p)[0] for q in self.extracted_ligands for p in q]
        bcvs = []

        for i, cav in enumerate(self.bcv):
            for j, bcv in enumerate(cav):
                try:
                    hr = HotspotReader(path=os.path.join(bcv, "out.zip")).read()
                    bcvs.append(hr)
                except:
                    bcvs.append(None)
                sg = Grid.get_single_grid(hr.super_grids, mask=False)

                for m in self.ligand_id:
                    for mol in m:
                        if mol in bcv:
                            for e in self.extracted_ligands:
                                for extracted in e:
                                    if mol in extracted:
                                        molecule = io.MoleculeReader(extracted)[0]
                                        molecule_grid = Grid.from_molecule(mol=molecule)

                                        perc_overlap = sg.percentage_overlap(molecule_grid)

                                        with open(self.all_overlaps[i+j], 'wb') as writer:
                                            writer.write(str(perc_overlap))

        # tasks

        bcv_by_score = {b.score(): b for b in [x for x in bcvs if x is not None]}
        print(bcv_by_score)
        top = bcv_by_score[sorted(bcv_by_score.keys(), reverse=True)[0]]
        sg = Grid.get_single_grid(top.super_grids, mask=False)

        for mol in mols:
            other = Grid.from_molecule(mol)
            po = sg.percentage_overlap(other=other)
            # output

            for path in self.overlap:
                if mol.identifier.split(":")[1][0:3] in path:
                    with open(path, "wb") as writer:
                        writer.write(str(po))

    def run(self, rerun=False):

        if not os.path.exists(self.biological_assembly) or rerun:
            self._download_pdb()
        try:
            if not os.path.exists(self.protonated) or rerun:
                self._protonate()
        except:
            pass


        if not os.path.exists(self.no_wat) or rerun:
            self._remove_wat_lig()

        for i, path in enumerate(self.other_pdbs):
            if not os.path.exists(path) or rerun:
                self._download_other_pdb(i)

        for i, path in enumerate(self.aligned_pdbs):
            paths = [os.path.exists(extracted)
                     for extracted in self.extracted_ligands[i]]
            paths.append(os.path.exists(path))

            if not all(paths) or rerun:
                self._align_other_to_main(i)

        for i, paths in enumerate(self.volumes):
            if not all([os.path.exists(p) for p in paths]) or rerun:
                self._get_ligand_volume(i)

        if not os.path.exists(self.buriedness) or rerun:
            self._get_buriedness_grid()

        self._get_cavities(min_vol=200)

        for i, path in enumerate(self.superstar):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                self._get_superstar(cav_id=i)

        for i, path in enumerate(self.hotspot):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                self._get_hotspot(cav_id=i)

        for i, paths in enumerate(self.bcv):
            if not all([os.path.exists(os.path.join(path, "out.zip")) for path in paths]) or rerun:
                self._get_bcv(cav_id=i)

        if not all([os.path.exists(overlap) for overlap in self.overlap] +
                   [os.path.exists(all_overlap) for all_overlap in self.all_overlaps]) or rerun:
            self._get_volume_overlap()


def main():

    # inputs
    df = pd.read_csv("readme.csv", index_col=0)
    buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']

    for index, row in df.iterrows():
        print row.apo, index
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

    hp = HotspotPipeline(pdb='1YES',
                         buriedness_method='ligsite',
                         protein_id=['2QFO', '2QG0'],
                         ligand_id=['A13_A51', 'A94'])

    hp.run(rerun=False)


if __name__ == '__main__':
    main()
