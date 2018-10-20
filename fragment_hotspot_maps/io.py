#
# This code is Copyright (C) 2015 The Cambridge Crystallographic Data Centre
# (CCDC) of 12 Union Road, Cambridge CB2 1EZ, UK and a proprietary work of CCDC.
# This code may not be used, reproduced, translated, modified, disassembled or
# copied, except in accordance with a valid licence agreement with CCDC and may
# not be disclosed or redistributed in any form, either in whole or in part, to
# any third party. All copies of this code made in accordance with a valid
# licence agreement as referred to above must contain this copyright notice.
#
# No representations, warranties, or liabilities are expressed or implied in the
# supply of this code by CCDC, its servants or agents, except where such
# exclusion or limitation is prohibited, void or unenforceable under governing
# law.
#
"""
The :mod:`fragment_hotspot_maps.io` module contains classes for the
reading and writing of hotspots.

The main classes of the :mod:`fragment_hotspot_maps.io` module are:

- :class:`fragment_hotspot_maps.io.HotspotWriter`
- :class:`fragment_hotspot_maps.io.HotspotReader`

TO DO:

"""
import zipfile
from os import listdir
from os.path import splitext, join, basename
import tempfile

from ccdc.protein import Protein
from ccdc import io

from grid_extension import Grid
from hotspot_calculation import HotspotResults


class PyMolWriter(object):
    """
    class to handle the writing of Hotspots to PyMol
    """
    def __init__(self, path):
        """

        :param path: directory (maybe zipped) containing hotspot information
        """
        self.path = path

    def write(self, hr, isosurface_threshold=[10, 14, 17], pharmacophores=false, island_labels=True):
        """hr result can be instance or list"""

class HotspotReader(object):
    """
    class to handle the reading of Hotspots
    """
    def __init__(self, path, sampled_probes=False):
        """

        :param path:
        """
        self._supported_interactions = ["apolar", "donor", "acceptor", "positive", "negative"]
        self._path = path

        ext = splitext(self._path)[1]
        if ext == ".zip":
            self._base = self._path_from_zip()
        else:
            self._base = path

        self._files = listdir(self._base)
        self._extensions = set([splitext(f)[1] for f in self._files if f != "" or f != ".py"])

        self.prot = Protein.from_file(join(self._base, [f for f in self._files if splitext(f)[1] == ".pdb"][0]))
        self.grid_dic, self.buriedness = self._get_grids()
        if sampled_probes:
            self.sampled_probes = self._get_sampled_probes()
        else:
            self.sampled_probes = None

        self.out_dir = None

    def _path_from_zip(self):
        """
        writes files to temp dir, returns base dir
        :param self:
        :return:
        """
        base = tempfile.mkdtemp()

        with zipfile.ZipFile(self._path) as hs_zip:
            hs_zip.extractall(base)
        d = splitext(basename(self._path))

        return join(base, d[0])

    def _get_grids(self):
        """
        create a grid dictorionary
        :return:
        """
        if ".grd" in self._extensions:
            grid_dic = {splitext(fname)[0]: Grid.from_file(join(self._base, fname))
                             for fname in [f for f in self._files
                                           if splitext(f)[1] == ".grd"
                                           and splitext(f)[0] in self._supported_interactions]}
            try:
                buriedness = Grid.from_file("buriedness.grd")
            except RuntimeError:
                buriedness = None

        elif ".dat" in self._extensions:
            grid_dic = {splitext(fname)[0]: Grid.from_array(join(self._base, fname))
                             for fname in [f for f in self._files
                                           if splitext(f)[1] == ".grd"
                                           and splitext(f)[0] in self._supported_interactions]}
            try:
                buriedness = Grid.from_array(join(self._base, "buriedness.dat"))
            except RuntimeError:
                buriedness = None

        else:
            raise IOError("grids not recognised")

        return grid_dic, buriedness

    def _get_sampled_probes(self):
        """
        retrieves sampled probes
        :return:
        """
        if ".mol2" not in self._extensions:
            return None
        else:

            try:
                probe_dict = {splitext(fname)[0].split("_")[0]: io.MoleculeReader(join(self._base,fname))
                             for fname in [f for f in self._files
                                           if splitext(f)[1] == ".mol2"
                                           and splitext(f)[0].split("_")[0] in self._supported_interactions]}
            except:
                 probe_dict = None

        return probe_dict

    def read(self):
        """
        creates a hotspot result
        :return:
        """
        return HotspotResults(protein=self.prot,
                              grid_dict=self.grid_dic,
                              fname=None,
                              sampled_probes=None,
                              buriedness=self.buriedness,
                              out_dir=self.out_dir)
    