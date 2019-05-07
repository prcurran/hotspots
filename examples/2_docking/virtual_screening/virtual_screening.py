"""
Example script for the Hotspot API manuscript
"""
from __future__ import print_function

import argparse
import os
import gzip
import shutil
import tempfile
from shutil import copyfile

from ccdc.docking import Docker
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.protein import Protein

from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.hs_docking import DockerSettings


class Organiser(argparse.ArgumentParser):
    """
    class to handle docking run with hotspot insights
    """

    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        # handle command line arguments
        self.add_argument(
            'path',
            help='path to directory containing DUD-E data.\n For example: "/home/usr/diverse/akt1"'
        )

        self.add_argument(
            'hotspot',
            help='path to directory containing Hotspot data.\n For example: "/home/usr/diverse/akt1"'
        )

        self.add_argument(
            'weight',
            help='penalty for not making the Hotspot protein hydrogen bond constraint"'
        )

        self.args = self.parse_args()
        self.hotspot = HotspotReader(os.path.join(self.args.hotspot, "hotspot", "out.zip")).read()
        for p, g in self.hotspot.super_grids.items():
            self.hotspot.super_grids[p] = g.max_value_of_neighbours()

    def dock(self):
        docker = Docker()
        docker.settings = DockerSettings()

        pfile = os.path.join(self.args.hotspot, "protein.mol2")
        with MoleculeWriter(pfile) as w:
            w.write(self.hotspot.protein)
        docker.settings.add_protein_file(pfile)

        docker.settings.binding_site = docker.settings.BindingSiteFromLigand(protein=docker.settings.proteins[0],
                                                                             ligand=MoleculeReader(
                                                                                 os.path.join(self.args.path,
                                                                                              "crystal_ligand.mol2"))[0])
        docker.settings.fitness_function = 'plp'
        docker.settings.autoscale = 1
        docker.settings.output_directory = tempfile.mkdtemp()
        docker.settings.output_file = "docked_ligands.mol2"

        with gzip.open(os.path.join(self.args.path, "actives_final.mol2.gz"), 'rb') as f_in:
            with open(os.path.join(self.args.path, "actives_final.mol2"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        with gzip.open(os.path.join(self.args.path, "decoys_final.mol2.gz"), 'rb') as f_in:
            with open(os.path.join(self.args.path, "decoys_final.mol2"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        docker.settings.add_ligand_file(os.path.join(self.args.path, "actives_final.mol2"), ndocks=5)
        docker.settings.add_ligand_file(os.path.join(self.args.path, "decoys_final.mol2"), ndocks=5)

        if int(self.args.weight) != 0:
            print(len(docker.settings.proteins[0].atoms))
            constraints = docker.settings.HotspotHBondConstraint.create(protein=docker.settings.proteins[0],
                                                                        hr=self.hotspot,
                                                                        weight=int(self.args.weight),
                                                                        min_hbond_score=0.05,
                                                                        max_constraints=1)

            print(constraints)
            for constraint in constraints:
                docker.settings.add_constraint(constraint)

        if not os.path.exists(os.path.join(self.args.path, "docking/{}".format(self.args.weight))):
            os.mkdir(os.path.join(self.args.path, "docking/{}".format(self.args.weight)))

        docker.dock(os.path.join(self.args.path, "docking/{}".format(self.args.weight), "hs_gold.conf"))
        results = Docker.Results(docker.settings)

        with MoleculeWriter(os.path.join(self.args.path, "docking/{}".format(self.args.weight), "docked_ligand.mol2")) as w:
            for d in results.ligands:
                w.write(d.molecule)

        copyfile(os.path.join(docker.settings.output_directory, "bestranking.lst"),
                 os.path.join(self.args.path, "docking/{}".format(self.args.weight), "bestranking.lst"))


def main():
    o = Organiser()
    o.dock()


if __name__ == "__main__":
    main()
