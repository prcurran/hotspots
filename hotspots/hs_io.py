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
The :mod:`hotspots.io` module contains classes for the
reading and writing of hotspots.

The main classes of the :mod:`hotspots.io` module are:

- :class:`hotspots.io.HotspotWriter`
- :class:`hotspots.io.HotspotReader`

TO DO:

"""
from __future__ import print_function
import zipfile
from os import listdir, walk
from os.path import splitext, join, basename, dirname, isdir
import tempfile
import shutil

from ccdc.molecule import Molecule, Atom
from ccdc.protein import Protein
from ccdc import io
import nglview as nv

from grid_extension import Grid
from calculation import Results
from template_strings import pymol_imports, pymol_arrow, pymol_protein, pymol_grids, pymol_display_settings, \
    pymol_load_zip, pymol_labels, pymol_mesh
from hs_utilities import Helper
from hs_pharmacophore import PharmacophoreModel


class HotspotWriter(Helper):
    """
    class to handle the writing of Hotspots to file. And creation of visualisation scripts.
    """
    class Settings(object):
        """class to hold writer settings"""
        def __init__(self):
            # format settings
            self.grid_extension = ".grd"

            # visual settings
            self.bg_color = "white"

            # protein
            self.surface = True
            self.surface_trim_factor = 13

            # ligand
            self.organic_sticks = True

            # grid
            self.isosurface_threshold=[10, 14, 17]
            self.grids = None
            self.transparency = 0.2
            self.grid_labels = True
            self.supported_grids = [".grd", ".ccp4", ".acnt"]

            #pharmacophore
            self.pharmacophore = False
            self.pharmacophore_labels = True
            self.pharmacophore_format = [".py"]

    def __init__(self, path, visualisation="pymol", grid_extension=".grd", zip_results=False, settings=None):
        """
        :param path: directory (maybe zipped) containing hotspot information
        """
        if settings is None:
            self.settings = self.Settings()
        else:
            self.settings = settings

        if grid_extension in self.settings.supported_grids:
            self.settings.grid_extension = grid_extension
        else:
            self.settings.grid_extension = ".grd"
            print("WARNING: Invalid grid file format provided. Default: '.grd' will be used")

        self.settings.visualisation = visualisation
        if visualisation == "ngl":
            self.settings.visualisation = visualisation
            if grid_extension != ".ccp4":
                self.settings.grid_extension = ".ccp4"
                print("WARNING: Grids must be .ccp4 format for visualisation with NGLViewer")

        self.path = path
        self.zipped = zip_results

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if traceback:
            print(traceback)

    def write(self, hr):
        """hr result can be instance or list"""

        if isinstance(hr, list):
            self.settings.grids = hr[0].super_grids.keys()
            self.container = "hotspot_boundaries"
            self.number_of_hotspots = len(hr)

            self.out_dir= Helper.get_out_dir(join(self.path, self.container))

            self._write_protein(hr[0].protein)
            if hr[0].pharmacophore:
                self.settings.pharmacophore = True
            #hts = [h.hotspot_result for h in hr]
            self._write_pymol(hr, self.zipped)

            for i, hotspot in enumerate(hr):
                self.out_dir = Helper.get_out_dir(join(self.path, self.container, str(i)))
                self.settings.isosurface_threshold = [round(hotspot.threshold,1)]

                bi = (Grid.super_grid(2, hotspot.best_island).max_value_of_neighbours()
                      > hotspot.threshold)

                self._write_grids(hotspot.super_grids, buriedness=None, mesh=bi)
                self._write_protein(hotspot.protein)

                if hotspot.pharmacophore:
                    self._write_pharmacophore(hotspot.pharmacophore)

                self._write_pymol(hotspot, False)

            self.out_dir = dirname(self.out_dir)
            if self.zipped:
                self.zip_results(join(dirname(self.out_dir), self.container))

        else:
            self.settings.grids = hr.super_grids.keys()
            self.container = "out"
            self.number_of_hotspots = 1

            self.out_dir = Helper.get_out_dir(join(self.path, self.container))
            self._write_grids(hr.super_grids, buriedness=hr.buriedness)
            self._write_protein(hr.protein)

            if hr.pharmacophore:
                self.settings.pharmacophore = True
                self._write_pharmacophore(hr.pharmacophore)
            self._write_pymol(hr, self.zipped)

            if self.zipped:
                self.zip_results(join(dirname(self.out_dir), self.container))

    def _write_grids(self, grid_dict, buriedness=None, mesh=None, out_dir=None):
        """
        writes grids to output directory
        :param grid_dict:
        :param buriedness:
        :return:
        """
        for p, g in grid_dict.items():
            fname = "{}{}".format(p, self.settings.grid_extension)
            if not out_dir:
                g.write(join(self.out_dir, fname))
            else:
                g.write(join(out_dir, fname))

        if buriedness:
            buriedness.write(join(self.out_dir, "buriedness{}".format(self.settings.grid_extension)))

        if mesh:
            mesh.write(join(self.out_dir, "mesh{}".format(self.settings.grid_extension)))

        if self.settings.grid_labels:
            labels = {"label_threshold_{}.mol2".format(threshold): self.get_label(grid_dict, threshold=threshold)
                      for threshold in self.settings.isosurface_threshold}

            for fname, label in labels.items():
                with io.MoleculeWriter(join(self.out_dir, fname)) as writer:
                    writer.write(label)

    def _write_protein(self, prot, out_dir=None):
        """
        writes protein to output directory
        :param prot:
        :return:
        """
        if not out_dir:
            out = join(self.out_dir, "protein.pdb")

        else:
            out = join(out_dir, "protein.pdb")

        with io.MoleculeWriter(out) as writer:
            writer.write(prot)

    def _write_pharmacophore(self, pharmacophore):
        """
        writes pharmacophore to output directory
        :param pharmacophore:
        :return:
        """
        out = [join(self.out_dir, "pharmacophore" + fmat) for fmat in self.settings.pharmacophore_format]
        for o in out:
            pharmacophore.write(o)

        label = self.get_label(pharmacophore)
        with io.MoleculeWriter(join(self.out_dir, "label_threshold_{}.mol2".format(pharmacophore.identifier))) \
                as writer:
            writer.write(label)

    def _write_pymol(self, hr, zipped=False):
        """
        Constructs PyMol python script to automatically display FH results
        :return:
        """
        pymol_out = pymol_imports()

        if self.settings.pharmacophore:
            pymol_out += pymol_arrow()

        if zipped:

            pymol_out += pymol_load_zip(basename(self.out_dir))

        pymol_out += pymol_protein(self.settings, self.zipped)

        if isinstance(hr, list):
            for i, h in enumerate(hr):
                self.settings.isosurface_threshold = [round(h.threshold, 1)]
                pymol_out += self._get_pymol_hotspot(h, i=i)
                pymol_out += pymol_mesh(i)

        else:
            pymol_out += self._get_pymol_hotspot(hr)

        pymol_out += pymol_display_settings(self.settings)

        if zipped:
            out = dirname(self.out_dir)
        else:
            out = self.out_dir

        with open('{}/pymol_results_file.py'.format(out), 'w') as pymol_file:
            pymol_file.write(pymol_out)

    def _get_pymol_hotspot(self, h, i=None):
        """
        core pymol sequence
        :param h: HotspotResult object
        :return:
        """
        pymol_out = ""
        if self.settings.grid_labels:

            for t in self.settings.isosurface_threshold:
                if i is not None:
                    fname = join(str(i), "label_threshold_{}.mol2".format(t))
                else:
                    fname = "label_threshold_{}.mol2".format(t)
                pymol_out += pymol_labels(fname= fname,
                                          objname="label_threshold_{}".format(t))

        pymol_out += pymol_grids(i, self.settings)

        if h.pharmacophore:
            if i is not None:
                f = join(str(i), "label_threshold_{}.mol2".format(h.pharmacophore.identifier))
            else:
                f = "label_threshold_{}.mol2".format(h.pharmacophore.identifier)

            pymol_out += h.pharmacophore._get_pymol_pharmacophore(lfile=f)

        return pymol_out

    def zip_results(self, archive_name, delete_directory=True):
        """
        Zips the output directory created for this :class:`hotspots.HotspotResults` instance, and
        removes the directory by default. The zipped file can be loaded directly into a new
        :class:`hotspots.HotspotResults` instance using the
        :func:`~hotspots.Hotspots.from_zip_dir` function

        :param archive_name: str, file path
        :param delete_directory: bool, remove the out directory once it has been zipped

        :return: None
        """
        self.archive_name = archive_name
        shutil.make_archive(self.archive_name, 'zip', self.out_dir)
        self.archive_loc = dirname("{}.zip".format(self.archive_name))
        if delete_directory:
            shutil.rmtree(self.out_dir)


class HotspotReader(object):
    """
    class to handle the reading of Hotspots
    """
    def __init__(self, path):
        """

        :param path:
        """
        self._supported_interactions = ["apolar", "donor", "acceptor", "positive", "negative"]
        self._supported_grids = [".grd", ".ccp4", ".acnt", ".dat"]
        self._not_hs_dir = ["best_islands", "peaks", "ins"]
        self._path = path

        ext = splitext(self._path)[1]
        if ext == ".zip":
            self._base = self._path_from_zip()
        else:
            self._base = path

        self._files = listdir(self._base)
        self._extensions = set([splitext(f)[1] for f in self._files if f != "" or f != ".py"])

        pfiles = [f for f in self._files if splitext(f)[1] == ".pdb"]

        if len(pfiles) > 1:
            print("WARNING! {} has been used as default protein".format(join(self._base, "protein.pdb")))
            pfiles = [p for p in self._files if f =="protein.pdb"]

        self.protein = Protein.from_file(join(self._base, pfiles[0]))
        self.hs_dir = [d for d in self._files
                       if isdir(join(self._base, d)) and d not in self._not_hs_dir]

    def _path_from_zip(self):
        """
        writes files to temp dir, returns base dir
        :param self:
        :return:
        """
        base = tempfile.mkdtemp()

        with zipfile.ZipFile(self._path) as hs_zip:
            hs_zip.extractall(base)
        #d = splitext(basename(self._path))

        return base

    def _get_grids(self, sub_dir=None):
        """
        create a grid dictorionary
        :return:
        """
        if sub_dir:
            base = join(self._base, sub_dir)
            self._files = listdir(base)
            self._extensions = set([splitext(f)[1] for f in self._files if f != '' or f != '.py'])
        else:
            base = self._base

        if ".dat" in self._extensions:
            grid_dic = {splitext(fname)[0]: Grid.from_array(join(base, fname))
                        for fname in [f for f in self._files
                                      if splitext(f)[1] == ".grd"
                                      and splitext(f)[0] in self._supported_interactions]}
            try:
                buriedness = Grid.from_array(join(self.base, "buriedness.dat"))
            except RuntimeError:
                buriedness = None

        else:
            ext = list(set(self._extensions).intersection(self._supported_grids))
            if len(ext) == 1:
                grid_dic = {splitext(fname)[0]: Grid.from_file(join(base, fname))
                                 for fname in [f for f in self._files
                                               if splitext(f)[1] == ext[0]
                                               and splitext(f)[0] in self._supported_interactions]}
                try:
                    buriedness = Grid.from_file("buriedness{}".format(ext[0]))
                except RuntimeError:
                    buriedness = None
            else:
                raise RuntimeError("Opps, something went wrong.")

        return grid_dic, buriedness

    def read(self, identifier=None):
        """
        creates a hotspot result
        :return:
        """

        if len(self.hs_dir) == 0:
            self.grid_dic, self.buriedness = self._get_grids()
            return Results(protein=self.protein,
                           super_grids=self.grid_dic,
                           buriedness=self.buriedness)

        else:
            hrs = []
            if identifier:
                self.grid_dic, self.buriedness = self._get_grids(sub_dir=str(identifier))
                return Results(protein=self.protein,
                               super_grids=self.grid_dic,
                               buriedness=self.buriedness)
            else:
                for dir in self.hs_dir:
                    self.grid_dic, self.buriedness = self._get_grids(sub_dir=dir)
                    hrs.append(Results(protein=self.protein,
                                       super_grids=self.grid_dic,
                                       buriedness=self.buriedness))
            return hrs





