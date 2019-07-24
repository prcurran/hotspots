import pickle
from urllib2 import urlopen as urlopen
from ccdc.protein import Protein
from ccdc.cavity import Cavity
from ccdc import io
import subprocess
from hotspots.grid_extension import Grid
from hotspots.atomic_hotspot_calculation import _AtomicHotspot, _AtomicHotspotResult
from hotspots.calculation import Runner, Buriedness, ExpBuriedness
from hotspots.result import Extractor, Results
from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.hs_utilities import Helper
import sys
import os
import shutil
import pandas as pd


def create_directory(path):
    """
    create a directory if it doesn't already exist

    :param path:
    :return: path
    """
    if not os.path.exists(path):
        os.mkdir(path)

    return path


def _align_proteins(reference, reference_chain, targets):
    """
    align proteins by chain

    :param `ccdc.protein.Protein` reference: align to me
    :param str reference_chain: align to this chain
    :param list targets: list of `ccdc.protein.Protein`
    :return tup: list(:class:`ccdc.protein.Protein`) and list (:classa`ccdc.molecule.Molecule`)
    """
    print("Aligning proteins to {}, chain {}...".format(reference.identifier, reference_chain))
    aligned_prots = []
    aligned_ligands = []

    reference = Protein.from_file(reference.fname)
    reference.add_hydrogens()

    for t in tqdm(targets):
        prot = Protein.from_file(t.fname)
        prot.detect_ligand_bonds()
        prot.add_hydrogens()
        chain = None
        for l in prot.ligands:
            if str(t.clustered_ligand) == str(l.identifier.split(":")[1][0:3]):
                try:
                    bs = Protein.BindingSiteFromMolecule(protein=prot,
                                                         molecule=l,
                                                         distance=6)
                    chain = bs.residues[0].identifier.split(":")[0]
                except:
                    chain = None

        if chain is None:
            print("\n        {} failed! No chain detected".format(t.identifier))
            continue
        try:
            binding_site_superposition = Protein.ChainSuperposition()
            (bs_rmsd, bs_transformation) = binding_site_superposition.superpose(reference[reference_chain],
                                                                                prot[chain])
            aligned_prots.append(prot)
            for lig in prot.ligands:
                if str(t.clustered_ligand) == str(lig.identifier.split(":")[1][0:3]):
                    if chain == str(lig.identifier.split(":")[0]):
                        aligned_ligands.append(lig)
        except:
            print("\n        {} failed!".format(t.identifier))
            continue

    return aligned_prots, aligned_ligands


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
        self.ligand_id = ligand_id

        # outputs

        # directories
        self.working_dir_base_base = create_directory("pdb_files/{}".format(self.pdb[1:3]))
        self.working_dir_base = create_directory("pdb_files/{}/{}".format(self.pdb[1:3], self.pdb))
        self.working_dir = create_directory("pdb_files/{}/{}/{}".format(self.pdb[1:3], self.pdb, self.buriedness_method))

        # files
        self.fragment_volumes = ["{}/volume_{}.pdb".format(self.working_dir_base, l) for l in ligand_id]


        self.biological_assembly = "{}/biological_assembly.pdb".format(self.working_dir_base)
        self.protonated = "{}/protonated.pdb".format(self.working_dir_base)
        self.no_wat = "{}/protonated_no_wat.pdb".format(self.working_dir_base)
        self.other_pdbs = ["{}/{}.pdb".format(self.working_dir_base, p) for p in self.protein_id]
        self.buriedness = "{}/buriedness_{}.grd".format(self.working_dir, self.buriedness_method)

        # these get set during the run depending on the number of cavities detected
        self.cavities = []
        self.superstar = []
        self.hotspot = []

    # def _get_ligand_volume(self):
    #     """
    #
    #
    #     :return:
    #     """
    #
    #     fragment_paths = ["/vagrant/github_pkg/hotspots/benchmark_set/{}/fragment.pdb".format(l)
    #                       for l in self.ligand_id]
    #
    #     # lead_paths = ["/vagrant/github_pkg/hotspots/benchmark_set/{}/lead.pdb".format(l)
    #     #               for l in self.ligand_id]
    #
    #     pass

    def _download_pdb(self):
        """
        download the pdb from the RCSB. pdb1 indicates the protein is the biologial assembly

        :return: None
        """

        # task
        url = "https://files.rcsb.org/download/{}.pdb1".format(self.pdb)
        print(url)
        pdbfile = urlopen(url).read()

        # output
        with open(self.biological_assembly, 'w') as out_file:
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
        pdbfile = urlopen(url).read()

        # output
        with open(self.other_pdbs[i], 'w') as out_file:
            out_file.write(pdbfile)

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

    # def _get_tract_map(self, cav_id, volume):
    #
    #     cav_path = os.path.join(self.working_dir, 'cavity_{}'.format(cav_id))
    #     cav_hotpot_path = os.path.join(cav_path, "hotspot", "out.zip")
    #     cav_tract_path = os.path.join(cav_path, "bcv_{}".format(volume), "out.zip")
    #
    #     if not os.path.exists(os.path.dirname(cav_tract_path)):
    #         os.mkdir(os.path.dirname(cav_tract_path))
    #
    #     hr = HotspotReader(cav_hotpot_path).read()
    #     e = Extractor(hr)
    #
    #     try:
    #         bcv = e.extract_volume(volume=volume)
    #         out_settings = HotspotWriter.Settings()
    #         out_settings.charged = False
    #         w = HotspotWriter(os.path.dirname(cav_tract_path), grid_extension=".grd", zip_results=True,
    #                           settings=out_settings)
    #
    #         w.write(bcv)
    #
    #     except Exception as e:
    #         with open(os.path.join(self.working_dir, 'cavity_{}'.format(cav_id),"bcv_{}".format(volume), 'failed.txt'), 'w') as w:
    #             w.write("Failed")

    # def run(self, rerun=False):
    #     if not os.path.exists(self.volume) or rerun:
    #         self._get_ligand_volume()

        if not os.path.exists(self.biological_assembly) or rerun:
            self._download_pdb()

        if not os.path.exists(self.protonated) or rerun:
            self._protonate()

        if not os.path.exists(self.no_wat) or rerun:
            self._remove_wat_lig()

        for i, pdb in enumerate(self.protein_id):
            if not os.path.exists(self.other_pdbs[i]) or rerun:
                self._download_other_pdb(i)


        # align to protein
        # extract ligand
        # calculate volume


        if not os.path.exists(self.buriedness) or rerun:
            self._get_buriedness_grid()

        self._get_cavities(min_vol=200)

        for i, path in enumerate(self.superstar):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                self._get_superstar(cav_id=i)

        for i, path in enumerate(self.hotspot):
            if not os.path.exists(os.path.join(path, "out.zip")) or rerun:
                self._get_hotspot(cav_id=i)


            #
            # for vol in volumes:
            #     if not os.path.exists(os.path.join(self.working_dir, "cavity_{}".format(i),
            #                                        "bcv_{}".format(vol), 'out.zip')) or rerun:
            #         self._get_tract_map(i, volume=vol)


def main():

    # inputs
    df = pd.read_csv("readme.csv", index_col=0)
    print(df)

    for index, row in df.iterrows():
        print row.apo, index
        hp = HotspotPipeline(pdb=row.apo,
                             buriedness_method='ligsite',
                             protein_id=[row.fragment, row.lead],
                             ligand_id=[row.fragment_ID, row.lead_ID]
        )
        hp.run(rerun=False)
        if index >= 1:
            break

if __name__ == '__main__':
    main()
