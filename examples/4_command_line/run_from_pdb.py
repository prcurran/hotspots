from __future__ import print_function

import argparse

from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter


class Organiser(argparse.ArgumentParser):
    """ class to initiate a Result instance from PDB file"""

    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        self.add_argument(
            'pdb',
            help='pdb code'
        )
        self.add_argument(
            'out_dir',
            help='output directory'
        )
        self.add_argument(
            '-b', '--buriedness_method',
            default='ghecom',
            help='method used to calculate buriedness (default = "ghecom")'
        )
        self.add_argument(
            '-z', '--zipped',
            default=True,
            help='compress output files(default = True)'
        )

        self.args = self.parse_args()

    def run(self):
        """from fragment hotspot calc from protein"""
        h = Runner()
        settings = Runner.Settings(sphere_maps=False)
        result = h.from_pdb(pdb_code=self.args.pdb,
                            charged_probes=True,
                            buriedness_method=self.args.buriedness_method,
                            nprocesses=5,
                            settings=settings)

        with HotspotWriter(path=self.args.out_dir, zip_results=self.args.zipped) as writer:
            writer.write(result)


def main():
    r = Organiser()
    r.run()


if __name__ == "__main__":
    main()
