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
        self.buriedness = os.path.join(self.working_dir, "buriedness_{}.grd")

        # these get set during the run depending on the number of cavities detected
        self.cavities = []
        self.superstar = []
        self.hotspot = []
        self.bcv = []

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
