"""
The :mod:`hotspots.hs_io` module was created to facilitate easy reading and
writing of Fragment Hotspot Map results.

There are multiple components to a :class:`hotspots.result.Result` including, the
protein, interaction grids and buriedness grid. It is therefore tedious to manually
read/write using the various class readers/writers. The Hotspots I/O organises this
for the user and can handle single :class:`hotspots.result.Result` or lists of
:class:`hotspots.result.Result`.


The main classes of the :mod:`hotspots.io` module are:

- :class:`hotspots.io.HotspotWriter`
- :class:`hotspots.io.HotspotReader`

"""
from __future__ import print_function

import shutil
import tempfile
import zipfile
from os import listdir
from os.path import splitext, join, basename, dirname, isdir

from ccdc import io
from ccdc.protein import Protein
from ccdc.molecule import Molecule, Atom
from hotspots.grid_extension import Grid
from hotspots.result import Results
from hotspots.hs_utilities import Helper
from hotspots.wrapper_pymol import PyMOLCommands, PyMOLFile


class HotspotWriter(Helper):
    """
    A class to handle the writing of a :class`hotspots.result.Result`. Additionally, creation of the
    PyMol visualisation scripts are handled here.

    :param str path: path to output directory
    :param str visualisation: "pymol" or "ngl" currently only PyMOL available
    :param str grid_extension: ".grd", ".ccp4" and ".acnt" supported
    :param bool zip_results: If True, the result directory will be compressed. (recommended)
    :param `hotspots.hs_io.HotspotWriter.Settings` settings: settings
    """

    class Settings(object):
        """
        A class to hold the :class:`hotspots.hs_io.HotspotWriter` settings
        """

        def __init__(self):
            self.grid_extension = ".grd"
            self.bg_color = "white"
            self.surface = False
            self.surface_trim_factor = 13
            self.organic_sticks = True
            self.isosurface_threshold = [10, 14, 17]
            self.grids = None
            self.transparency = 0.2
            self.grid_labels = True
            self.supported_grids = [".grd", ".ccp4", ".acnt"]
            self.pharmacophore = False
            self.pharmacophore_labels = True
            self.pharmacophore_format = [".py"]
            self.container = 'out'
            self.identifier_tag = [0]
            self.colour_dict = {'acceptor':'red',
                                 'donor':'blue',
                                 'apolar':'yellow',
                                 'negative':'purple',
                                 'positive':'cyan'}

    def __init__(self, path, grid_extension=".grd", zip_results=True, settings=None):
        if settings is None:
            self.settings = self.Settings()
        else:
            self.settings = settings

        if grid_extension in self.settings.supported_grids:
            self.settings.grid_extension = grid_extension
        else:
            self.settings.grid_extension = ".grd"
            print("WARNING: Invalid grid file format provided. Default: '.grd' will be used")

        self.path = self.get_out_dir(path)

        self.zip_results = zip_results

        self.pymol_out = PyMOLFile()
        if self.zip_results:
            # unzip and change working directory
            self.pymol_out.commands += PyMOLCommands.unzip_dir(self.settings.container)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if traceback:
            print(traceback)

    def write(self, hr):
        """
        writes the Fragment Hotspot Maps result to the output directory and create the pymol visualisation file

        :param hr: a Fragment Hotspot Maps result or list of results
        :type hr: `hotspots.result.Result`

        >>> from hotspots.calculation import Runner
        >>> from hotspots.hs_io import HotspotWriter

        >>> r = Runner
        >>> result = r.from_pdb("1hcl")
        >>> out_dir = <path_to_out>
        >>> with HotspotWriter(out_dir) as w:
        >>>     w.write(result)

        """
        container = Helper.get_out_dir(join(self.path, self.settings.container))

        if isinstance(hr, list):
            if len({h.identifier for h in hr}) != len(hr):
                # if there are not unique identifiers, create some.
                for i, h in enumerate(hr):
                    h.identifier = f"hotspot_{i}"
            for h in hr:
                self._single_write(container, h)

        else:
            if not hr.identifier:
                hr.identifier = "hotspot"
            self._single_write(container, hr)

        self._write_pymol_isoslider(hr)
        self.pymol_out.commands += PyMOLCommands.display_settings(self.settings)
        self.pymol_out.commands += PyMOLCommands.push_to_wd()

        if self.zip_results:
            self.compress()

        self.pymol_out.write(join(self.path, "pymol_file.py"))

    def _single_write(self, path, hr):
        hr.out_dir = Helper.get_out_dir(join(path, hr.identifier))
        self._write_grids(hr.out_dir, hr.super_grids, buriedness=hr.buriedness)
        self._write_protein(hr.out_dir, hr.protein)

        relpath = f'{hr.identifier}'
        self._write_pymol_objects(relpath, hr)

    def _write_grids(self, path, grid_dict, buriedness=None):
        """
        Write probe grids to output directory

        :param grid_dict: hotspot result grid dict
        :param buriedness: buriedness grid

        :type grid_dict: hotspot.grid_extension.Grid
        :type buriedness: hotspot.grid_extension.Grid
        """
        for p, g in grid_dict.items():
            g.write(join(path, f"{p}{self.settings.grid_extension}"))

        if buriedness:
            buriedness.write(join(path, f"buriedness{self.settings.grid_extension}"))

    @staticmethod
    def _write_protein(path, prot):
        """
        writes protein to output directory

        :param prot: a protein
        :type prot: `ccdc.protein.Protein`

        """
        with io.MoleculeWriter(join(path, "protein.pdb")) as writer:
            writer.write(prot)

    def _write_pymol_isoslider(self, hr):
        """
        generates the commands for an isoslider

        :param hr: a hotspot result
        :type hr: `hotspots.results.Results`
        """
        # Isosurface obj's take the name: "surface_{hotspots ID}_{probe ID}"
        # e.g. "surface_hotspotA_apolar"
        if isinstance(hr, list):
            surface_dic = {h.identifier: [f"surface_{g}_{h.identifier}"
                                          for g in h.super_grids.keys()] for h in hr}
            surface_value_dic = {h.identifier: [g for g in h.super_grids.values()] for h in hr}
        else:
            surface_dic = {hr.identifier: [f"surface_{g}_{hr.identifier}"
                                           for g in hr.super_grids.keys()]}
            surface_value_dic = {hr.identifier: [g for g in hr.super_grids.values()]}

        # surface_dict.values() = [['g0', 'g1',...], ['g0',...]]
        max_value = max([round(g.extrema[1], 1) for lst in surface_value_dic.values() for g in lst])
        min_value = 0
        self.pymol_out.commands += PyMOLCommands.isoslider(surface_dic, min_value, max_value)

    def _write_pymol_objects(self, relpath, hr, load_prot=True):
        """
        generates pymol commands associated with an indivdual hotspot

        :param relpath: path to the directory holding associated files
        :param hr: hotspot result

        :type relpath: str
        :type hr: `hotspots.results.Results`

        """
        default_level = 5
        # load grids and create isosurfaces
        for p in hr.super_grids.keys():
            objname = f'{p}_{hr.identifier}'
            self.pymol_out.commands += PyMOLCommands.load(fname=f'{relpath}/{p}{self.settings.grid_extension}',
                                                          objname=objname)

            # surface_10_apolar_hotspotid
            surface_objname = f'surface_{objname}'
            self.pymol_out.commands += PyMOLCommands.isosurface(grd_name=objname,
                                                                isosurface_name=surface_objname,
                                                                level=default_level,
                                                                color=self.settings.colour_dict[p])

            self.pymol_out.commands += PyMOLCommands.pymol_set(setting_name='transparency',
                                                               value=self.settings.transparency,
                                                               selection=surface_objname)

        group_members = [f'{p}_{hr.identifier}' for p in hr.super_grids.keys()] + \
                        [f'surface_{p}_{hr.identifier}' for p in hr.super_grids.keys()]

        self.pymol_out.commands += PyMOLCommands.group(group_name=hr.identifier,
                                                       members=group_members)

        # generate grid labels
        labels = hr.grid_labels()

        for p, dic in labels.items():
            i = 0
            group_me = []
            for coord, value in dic.items():
                objname = f"PS_{p}_{i}"
                group_me.append(objname)
                self.pymol_out.commands += PyMOLCommands.pseudoatom(objname=objname,
                                                                    coords=coord,
                                                                    label=f'{round(value, 1)}')
                group_me.append(objname)
                i += 1
            self.pymol_out.commands += PyMOLCommands.group(f'label_{p}_{hr.identifier}', group_me)

        self.pymol_out.commands += PyMOLCommands.group("labels", [f'label_{p}_{hr.identifier}'
                                                                  for p in hr.super_grids.keys()])

        # load up the protein
        if load_prot:
            self.pymol_out.commands += PyMOLCommands.load(f'{relpath}/protein.pdb', f'protein_{hr.identifier}')

        # find contributing residues

    def compress(self, delete_directory: bool = True) -> None:
        """
        compresses the output directory created for this :class:`hotspots.HotspotResults` instance, and
        removes the directory by default. The zipped file can be loaded directly into a new
        :class:`hotspots.HotspotResults` instance using the
        :func:`~hotspots.Hotspots.from_zip_dir` function

        :param bool delete_directory: remove the out directory once it has been zipped
        """
        out_dir = join(self.path, self.settings.container)
        shutil.make_archive(out_dir, 'zip', out_dir)
        if delete_directory:
            shutil.rmtree(join(self.path, self.settings.container))


class HotspotReader(object):
    """
    A class to organise the reading of a :class:`hotspots.result.Result`

    :param str path: path to the result directory (can be .zip directory)
    """

    def __init__(self, path):
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
            pfiles = [p for p in self._files if f == "protein.pdb"]

        self.protein = Protein.from_file(join(self._base, pfiles[0]))
        self.hs_dir = [d for d in self._files
                       if isdir(join(self._base, d)) and d not in self._not_hs_dir]

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        try:
            shutil.rmtree(self._base)
        except:
            pass

    def _path_from_zip(self):
        """
        writes files to temp dir, returns base dir
        :param self:
        :return:
        """
        base = tempfile.mkdtemp()

        with zipfile.ZipFile(self._path) as hs_zip:
            hs_zip.extractall(base)
        # d = splitext(basename(self._path))

        return base

    def _get_grids(self, sub_dir=None):
        """
        create a grid dic
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
        creates a single or list of :class:`hotspots.result.Result` instance(s)

        :param str identifier: for directories containing multiple Fragment Hotspot Map results,
        identifier is the subdirectory for which a :class:`hotspots.result.Result` is requried

        :return: `hotspots.result.Result` a Fragment Hotspot Map result

        >>> from hotspots.hs_io import HotspotReader

        >>> path = "<path_to_results_directory>"
        >>> result = HotspotReader(path).read()


        """
        if len(self.hs_dir) == 0:
            self.grid_dic, self.buriedness = self._get_grids()
            shutil.rmtree(self._base)
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

            shutil.rmtree(self._base)
            return hrs
