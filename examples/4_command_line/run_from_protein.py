from __future__ import print_function

import argparse
from os.path import dirname, abspath

from ccdc.cavity import Cavity
from ccdc.protein import Protein

from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter


class Organiser(argparse.ArgumentParser):
    """ class to initiate a HotspotResults instance"""

    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        self.add_argument(
            'prot_fname',
            help='pdb file path'
        )
        self.add_argument(
            '-b', '--buriedness_method',
            default='ghecom',
            help='method used to calculate buriedness (default = "ghecom")'
        )

        self.add_argument(
            '-w', '--remove_waters',
            default=True,
            help='remove waters in prot preparation(default = True)'
        )
        
        self.add_argument(
            '-p', '--prepare',
            default=True,
            help='prepare_protein(default = True)'
        )

        self.add_argument(
            '-z', '--zipped',
            default=True,
            help='compress output files(default = True)'
        )

        self.args = self.parse_args()

        self.in_dir = dirname(abspath(self.args.prot_fname))

    def prepare_protein(self):
        """default protein preparation settings on the protein"""
        self.prot = Protein.from_file(self.args.prot_fname)
        self.prot.add_hydrogens()
        for lig in self.prot.ligands:
            self.prot.remove_ligand(lig.identifier)
        self.prot.remove_all_metals()
        if self.args.remove_waters:
            self.prot.remove_all_waters()

    def run(self, cavity=True):
        """from fragment hotspot calc from protein"""
        h = Runner()
        settings = Runner.Settings(sphere_maps=False)
        if self.args.prepare:
            self.prepare_protein()
        else:
            self.prot = Protein.from_file(self.args.prot_fname)

        if cavity:
            cavs = Cavity.from_pdb_file(self.args.prot_fname)
            print(cavs)
        else:
            cavs = None

        result = h.from_protein(protein=self.prot,
                                charged_probes=False,
                                buriedness_method=self.args.buriedness_method,
                                cavities=cavs,
                                nprocesses=5,
                                settings=settings)

        with HotspotWriter(path=self.in_dir, zip_results=self.args.zipped) as writer:
            writer.write(result)


def main():
    r = Organiser()
    r.run(cavity=False)


if __name__ == "__main__":
    main()
