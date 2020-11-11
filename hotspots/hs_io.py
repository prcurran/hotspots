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
from os import listdir, walk, getenv
from os.path import splitext, join, isdir, isfile, basename

from ccdc import io
from ccdc.protein import Protein
from ccdc.utilities import PushDir
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
            self.output_superstar = False
            self.output_weighted = False
            self.output_buriedness = True
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
            # e.g. self.protein_color_dic = {f"protein_{hr.identifier}": "slate"}
            self.protein_color_dic = {}
            self.colour_dict = {'acceptor':'red',
                                 'donor':'blue',
                                 'apolar':'yellow',
                                 'negative':'purple',
                                 'positive':'cyan',
                                    'buriedness': 'gray'}

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
            print(hr)
            if len({h.identifier for h in hr}) != len(hr):
                # if there are not unique identifiers, create some.
                for i, h in enumerate(hr):
                    h.identifier = f"hotspot-{i}"
            for h in hr:
                self._single_write(container, h)

        else:
            if not hr.identifier:
                hr.identifier = "hotspot"
            self._single_write(container, hr)

        self._write_pymol_isoslider(hr)
        self.pymol_out.commands += PyMOLCommands.background_color(self.settings.bg_color)
        self.pymol_out.commands += PyMOLCommands.push_to_wd()

        if self.zip_results:
            self.compress()

        self.pymol_out.write(join(self.path, "pymol_file.py"))

    def _single_write(self, path, hr):
        hr.out_dir = Helper.get_out_dir(join(path, hr.identifier))

        self._write_grids(hr)
        self._write_protein(hr.out_dir, hr.protein)

        relpath = f'{hr.identifier}'
        self._write_pymol_objects(relpath, hr)

    def _write_grids(self, hr):
        """
        Write probe grids to output directory

        :param grid_dict: hotspot result grid dict
        :param buriedness: buriedness grid

        :type grid_dict: hotspot.grid_extension.Grid
        :type buriedness: hotspot.grid_extension.Grid
        """
        for p, g in hr.super_grids.items():
            g.write(join(hr.out_dir, f"{p}{self.settings.grid_extension}"))

        if self.settings.output_buriedness and hr.buriedness:
            hr.buriedness.write(join(hr.out_dir, f"buriedness{self.settings.grid_extension}"))

        if self.settings.output_superstar and hr.superstar:
            for p, g in hr.superstar.items():
                g.write(join(hr.out_dir, f"superstar_{p}{self.settings.grid_extension}"))

        if self.settings.output_weighted and hr.weighted_superstar:
            for p, g in hr.weighted_superstar.items():
                g.write(join(hr.out_dir, f"weighted_{p}{self.settings.grid_extension}"))

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
        # Isosurface obj's take the name: "surface_{probe ID}_{hotpsot ID}"
        # e.g. "surface_apolar_hotspotA"
        if not isinstance(hr, list):
            hr = [hr]

        # the hotspot grids are always output
        surface_dic = {h.identifier: {'fhm': [f"surface_{g}_{h.identifier}" for g in h.super_grids.keys()]}
                       for h in hr}

        surface_value_dic = {h.identifier: {"fhm": max([round(g.extrema[1], 1) for g in h.super_grids.values()])}
                             for h in hr}

        for h in hr:
            if self.settings.output_superstar and h.superstar:
                surface_dic[h.identifier].update({'superstar': [f"surface_superstar_{g}_{h.identifier}"
                                                                for g in h.superstar.keys()]})

                surface_value_dic[h.identifier].update({'superstar': max([round(g.extrema[1], 1)
                                                                          for g in h.superstar.values()])})

            if self.settings.output_weighted and h.weighted_superstar:
                surface_dic[h.identifier].update({'weighted': [f"surface_weighted_superstar_{g}_{h.identifier}"
                                                                for g in h.weighted_superstar.keys()]})

                surface_value_dic[h.identifier].update({'weighted': max([round(g.extrema[1], 1)
                                                                          for g in h.weighted_superstar.values()])})

            if self.settings.output_buriedness and h.buriedness:
                surface_dic[h.identifier].update({'buriedness': [f"surface_buriedness_{h.identifier}"]})

                surface_value_dic[h.identifier].update({'buriedness': 8})

        min_value = 0
        print(surface_dic)
        print(surface_value_dic)

        self.pymol_out.commands += PyMOLCommands.isoslider(surface_dic, surface_value_dic)

    def _write_pymol_isosurfaces(self, dict, relpath, identifier, dict_type):
        """
        Loads grids and generates isosurfaces

        :param dict: interaction grid dictionary
        :param relpath: result containing directory
        :param identifier: hotspot identifier
        :param dict_type: superstar, fhm or weighted_superstar

        :type dict: dict
        :type relpath: str
        :type identifier: str
        :type dict_type: str

        :return: pymol commands
        :rtype: str
        """
        cmd = ""
        default_level = 5
        # load grids and create isosurfaces
        group_members = []
        for p in dict.keys():
            if dict_type == 'fhm':
                objname = f'{p}_{identifier}'
                fname = f'{relpath}/{p}{self.settings.grid_extension}'
            else:
                objname = f'{dict_type}_{p}_{identifier}'
                fname = f'{relpath}/{dict_type}_{p}{self.settings.grid_extension}'
            cmd += PyMOLCommands.load(fname=fname, objname=objname)

            # surface_10_apolar_hotspotid
            surface_objname = f'surface_{objname}'
            cmd += PyMOLCommands.isosurface(grd_name=objname,
                                            isosurface_name=surface_objname,
                                            level=default_level,
                                            color=self.settings.colour_dict[p])

            cmd += PyMOLCommands.pymol_set(setting_name='transparency',
                                           value=self.settings.transparency,
                                           selection=surface_objname)

            group_members.extend([objname, f"surface_{objname}"])

        cmd += PyMOLCommands.group(group_name=identifier, members=group_members)
        return cmd

    def _write_pymol_objects(self, relpath, hr, load_prot=True):
        """
        generates pymol commands associated with an indivdual hotspot

        :param relpath: path to the directory holding associated files
        :param hr: hotspot result

        :type relpath: str
        :type hr: `hotspots.results.Results`

        """
        self.pymol_out.commands += self._write_pymol_isosurfaces(hr.super_grids, relpath, hr.identifier, "fhm")

        if self.settings.output_superstar and hr.superstar:
            self.pymol_out.commands += self._write_pymol_isosurfaces(hr.superstar, relpath, hr.identifier, "superstar")

        if self.settings.output_weighted and hr.weighted_superstar:
            self.pymol_out.commands += self._write_pymol_isosurfaces(hr.weighted_superstar, relpath, hr.identifier, "weighted")

        if self.settings.output_buriedness and hr.buriedness:
            default_level=3
            objname = f'buriedness_{hr.identifier}'
            fname = f'{relpath}/buriedness{self.settings.grid_extension}'

            self.pymol_out.commands += PyMOLCommands.load(fname=fname, objname=objname)

            # surface_10_apolar_hotspotid
            surface_objname = f'surface_{objname}'
            self.pymol_out.commands += PyMOLCommands.isosurface(grd_name=objname,
                                                       isosurface_name=surface_objname,
                                                       level=default_level,
                                                       color=self.settings.colour_dict["buriedness"])

            self.pymol_out.commands += PyMOLCommands.pymol_set(setting_name='transparency',
                                                      value=self.settings.transparency,
                                                      selection=surface_objname)

        group_members = [f'buriedness_{hr.identifier}', f'surface_buriedness_{hr.identifier}']

        self.pymol_out.commands += PyMOLCommands.group(group_name=hr.identifier, members=group_members)

        # generate grid labels
        labels = hr.grid_labels()

        for p, dic in labels.items():
            i = 0
            group_me = []
            for coord, value in dic.items():
                objname = f"PS_{p}_{hr.identifier}_{i}"
                group_me.append(objname)
                self.pymol_out.commands += PyMOLCommands.pseudoatom(objname=objname,
                                                                    coords=coord,
                                                                    label=f'{round(value, 1)}')
                group_me.append(objname)
                i += 1
            self.pymol_out.commands += PyMOLCommands.group(f'label_{p}_{hr.identifier}', group_me)

        self.pymol_out.commands += PyMOLCommands.group(f"labels_{hr.identifier}", [f'label_{p}_{hr.identifier}'
                                                                  for p in hr.super_grids.keys()])

        # load up the protein
        if load_prot:
            self.pymol_out.commands += PyMOLCommands.load(f'{relpath}/protein.pdb', f'protein_{hr.identifier}')
            if len(self.settings.protein_color_dic) > 0:
                self.pymol_out += PyMOLCommands.color("slate", f'protein_{hr.identifier}')
            self.pymol_out.commands += PyMOLCommands.show("cartoon", f'protein_{hr.identifier}')
            self.pymol_out.commands += PyMOLCommands.hide("line", f'protein_{hr.identifier}')
            self.pymol_out.commands += PyMOLCommands.show("sticks", "organic")

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
    def __init__(self, path, read_superstar=False, read_weighted=False, read_buriedness=True):
        self.read_superstar = read_superstar
        self.read_weighted = read_weighted
        self.read_buriedness = read_buriedness
        self.supported_interactions = {"apolar", "donor", "acceptor", "positive", "negative"}
        self.supported_grid_extensions = {"grd", "ccp4", "acnt", "dat"}
        self.supported_protein_extensions = "pdb"
        self.path = path

    def __enter__(self):
        self.top_extension = splitext(self.path)[1]

        if self.top_extension == ".zip":
            temp = tempfile.mkdtemp(prefix=getenv('TMPDIR_PREFIX', None))
            print("extract to ", temp)
            with zipfile.ZipFile(self.path) as hs_zip:
                hs_zip.extractall(temp)

            self.path = temp

        else:
            self.path = self.path

        return self

    def __exit__(self, type, value, traceback):
        if self.top_extension == ".zip":
            try:
                shutil.rmtree(self.path)
            except:
                pass

    def _generate_result(self, path):
        with PushDir(path):
            files = set(listdir(path))

            # fetch protein - this should always be protein.pdb
            prot_name = [f for f in files if f.split(".")[1] == self.supported_protein_extensions][0]
            prot = Protein.from_file(prot_name)
            files.remove(prot_name)

            # there should only be one grid extension in the directory, if there are more
            # then you can manually read in your results
            grid_extension = {f.split(".")[1] for f in files}.intersection(self.supported_grid_extensions)
            if len(grid_extension) > 1:
                raise IndexError("Too many grid types, create `hotspots.result.Results` manually")

            elif len(grid_extension) < 1:
                raise IndexError("No supported grid types found")

            elif list(grid_extension)[0] == "dat":
                raise NotImplementedError("Will put this in if requested")

            else:
                grid_extension = list(grid_extension)[0]

            # read hotspot grids
            stripped_files = {f.split(".")[0] for f in files}
            hotspot_grids = stripped_files.intersection(self.supported_interactions)
            super_grids = {p: Grid.from_file(f"{p}.{grid_extension}") for p in hotspot_grids}

            # read superstar grids
            if len([f.startswith("superstar") for f in files]) > 0 and self.read_superstar:
                superstar_grids = {p: Grid.from_file(f"superstar_{p}.{grid_extension}") for p in hotspot_grids}
            else:
                superstar_grids = None

            # read weighted_superstar grids
            if len([f.startswith("weighted") for f in files]) > 0 and self.read_weighted:
                weighted_grids = {p: Grid.from_file(f"weighted_{p}.{grid_extension}") for p in hotspot_grids}
            else:
                weighted_grids = None

            # fetch buriedness grid
            try:
                buriedness_name = [f for f in files if f.startswith("buriedness")][0]
            except IndexError:
                buriedness_name = None

            if buriedness_name and self.read_buriedness:
                buriedness = Grid.from_file(buriedness_name)
            else:
                buriedness = None

        return Results(super_grids=super_grids,
                       protein=prot,
                       buriedness=buriedness,
                       superstar=superstar_grids,
                       weighted_superstar=weighted_grids,
                       identifier=basename(path))

    def read(self, identifier=None):
        """
        creates a single or list of :class:`hotspots.result.Result` instance(s)

        :param str identifier: for directories containing multiple Fragment Hotspot Map results,
        identifier is the subdirectory for which a :class:`hotspots.result.Result` is requried

        :return: `hotspots.result.Result` a Fragment Hotspot Map result

        >>> from hotspots.hs_io import HotspotReader

        >>> path = "<path_to_results_directory>"
        >>> with HotspotReader(path) as reader:
        >>>    reader.read()

        """
        root, d, files = list(walk(self.path))[0]

        if len(d) == 0:
            # old style single result all files in os.listdir(self.base)
            hr = self._generate_result(path=self.path)

        elif len(d) == 1:
            # new style single result
            hr = self._generate_result(path=join(self.path, d[0]))

        else:
            # more than one hotspot
            if identifier:
                hr = self._generate_result(path=join(self.path, identifier))
            else:
                hr = [self._generate_result(path=join(self.path, x)) for x in d]

        return hr
