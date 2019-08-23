"""
Atomic Hotspot detection is the first step in the Fragment Hotspot Maps algorithm and is implemented through SuperStar.

The main class of the :mod:`hotspots.atomic_hotspot_calculation` module are:

- :class:`hotspots.atomic_hotspot_calculation.AtomicHotspot`

More information about the SuperStar method is available:
    - SuperStar: A knowledge-based approach for identifying interaction sites in proteins M. L. Verdonk, J. C. Cole and R. Taylor, J. Mol. Biol., 289, 1093-1108, 1999 [DOI: 10.1006/jmbi.1999.2809]
"""
from __future__ import print_function

import glob
import subprocess
import sys
import tempfile
from os import environ, mkdir
from os.path import join, exists, isfile, dirname

from ccdc.io import csd_directory, MoleculeWriter
from ccdc.utilities import PushDir
from concurrent import futures

try:
    from grid_extension import Grid
except ImportError:
    from ccdc.utilities import Grid


def _run_job(args):
    """
    private method

    Runs an Atomic Hotspot calculation using the SuperStar algorithm.
    :param tup or list args: Command, Jobname, SuperStar Environment Variable and Temporary directory
    :return: tup, temporary directory and jobname
    """
    cmd, jobname, superstar_env, temp_dir = args
    print(temp_dir)
    env = environ.copy()
    env.update(superstar_env)
    with PushDir(temp_dir):
        subprocess.call(cmd, shell=sys.platform != 'win32', env=env)
    return temp_dir, jobname


class _AtomicHotspot(object):
    """
    A class for handling the calculation of Atomic Hotspots using SuperStar
    """
    class Settings(object):
        """ handles the adjustable settings of Atomic Hotspot calculation.

        NB: not all settings are exposed here, for other settings please look at
        `hotspots.atomic_hotspot_calculation._atomic_hotspot_ins()`

        :param database: database from which the underlying interaction libraries are extracted
        :param map_background_value: the default value on the output grid
        :param box_border: padding on the output grid
        :param min_propensity: the minimum propensity value that can be assigned to a grid. value = 1 is random
        :param superstar_sigma: the sigma value for the gaussian smoothing function
        """

        def __init__(self, database='CSD', map_background_value=1, box_border=10, min_propensity=0,
                     superstar_sigma=0.5):
            self.database = database
            self.mapbackgroundvalue = map_background_value
            self.boxborder = box_border
            self.minpropensity = min_propensity
            self.superstar_sigma = superstar_sigma
            self.superstar_executable, self.superstar_env = self._set_environment_variables()
            self.temp_dir = tempfile.mkdtemp()
            self._csd_atomic_probes = {}
            self._pdb_atomic_probes = {}

        @staticmethod
        def _set_environment_variables():
            """
            private method

            sets up superstar environment variables
            :return: superstar executable (str), superstar env(str)
            """
            base = csd_directory()
            main_dir = environ.get('MAINDIR')
            if main_dir:
                if sys.platform == 'win32':
                    superstar_executable = 'superstar_app.exe'
                else:
                    superstar_executable = ' '.join([join(environ['MAINDIR'], 'run.sh'), 'superstar_app.x'])
                superstar_env = dict()
            else:
                if sys.platform == 'win32':
                    base = dirname(base)
                    merc = glob.glob(join(base, 'mercury*'))
                    if type(merc) is list:
                        try:
                            merc = merc[0]
                        except IndexError:
                            raise IndexError("No mercury path found, check API version")

                    superstar_executable = join(merc, 'superstar_app.exe')
                    if not isfile(superstar_executable):
                        superstar_executable = join(merc, 'superstar.exe')
                        if not isfile(superstar_executable):
                            raise IOError("superstar executable not found")

                    superstar_env = dict(SUPERSTAR_ISODIR=str(join(base, 'isostar_files', 'istr')),
                                         SUPERSTAR_ROOT=str(join(base, "Mercury"))
                                         )

                elif sys.platform == 'darwin':
                    print("OS X not supported")

                else:
                    base = dirname(base)
                    superstar_executable = join(base, 'bin', 'superstar')
                    superstar_env = dict(SUPERSTAR_ISODIR=str(join(base, 'isostar_files', 'istr')),
                                         SUPERSTAR_ROOT=str(base)
                                         )

            return superstar_executable, superstar_env

        @property
        def atomic_probes(self):
            """
            The probe atoms to be used by the Atomic Hotspot calculation.

            Available atomic probes:
                for the CSD:
                    - "Alcohol Oxygen",
                    - "Water Oxygen",
                    - "Carbonyl Oxygen",
                    - "Oxygen Atom",
                    - "Uncharged NH Nitrogen",
                    - "Charged NH Nitrogen",
                    - "RNH3 Nitrogen",
                    - "Methyl Carbon",
                    - "Aromatic CH Carbon",
                    - "C-Cl Chlorine"
                    - "C-F Fluorine",
                    - "Cyano Nitrogen",
                    - "Sulphur Atom",
                    - "Chloride Anion",
                    - "Iodide Anion"

                for the PDB:
                    - "Alcohol Oxygen",
                    - "Water Oxygen",
                    - "Carbonyl Oxygen",
                    - "Amino Nitrogen",
                    - "Aliphatic CH Carbon"
                    - "Aromatic CH Carbon"

            :return dict: key= "identifier", value="SuperStar identifier"

            """
            if self.database == 'CSD':
                return self._csd_atomic_probes

            elif self.database == 'PDB':
                return self._pdb_atomic_probes

            else:
                raise TypeError("Database must be 'CSD' or 'PDB'")

        @atomic_probes.setter
        def atomic_probes(self, probes):
            if self.database == 'CSD':
                self._csd_atomic_probes.update(probes)

            elif self.database == 'PDB':
                self._pdb_atomic_probes.update(probes)

            else:
                raise TypeError("Database must be 'CSD' or 'PDB'")

    class InstructionFile(object):
        """ handles the instruction files for the commandline call of the Atomic Hotspot calculation

        :param str jobname: general format "<probename>.ins"
        :param str probename: identifier for the atomic probe
        :param `hotspots.AtomicHotspot.Settings` settings: supplied if settings are adjusted from default
        :param tup cavity: (float(x), float(y), float(z)) describing the centre of a cavity
        """

        def __init__(self, jobname, probename, settings, cavity=None):
            self.ins_str = _atomic_hotspot_ins(jobname, probename, settings)
            if cavity:
                self.ins_str += '\nCAVITY_ORIGIN {} {} {}'.format(round(cavity[0], 1),
                                                                  round(cavity[1], 1),
                                                                  round(cavity[2], 1)
                                                                  )
            else:
                self.ins_str += '\nSUBSTRUCTURE ALL'

        def write(self, fname):
            """
            writes out the instruction file, enables calculation to be repeated

            :param str fname: path of the output file
            """
            with open(fname, "w") as f:
                f.write(self.ins_str)

    def __init__(self, settings=None):
        if not settings:
            self.settings = self.Settings()

        else:
            self.settings = settings

    def _get_cmd(self, protein, cavity_origin, out=None):
        """
        private method

        constructs the commandline str required to run atomic hotspot calculation
        :param str jobname: general format (<probe_type>.ins)
        :param str probename: probe identifier
        :param str out: output directory (outputting ins maybe useful for debugging)
        :return:
        """
        cmds = []
        if not out:
            out = self.settings.temp_dir

        with PushDir(out):
            with MoleculeWriter(join(out, 'protein.pdb')) as writer:
                writer.write(protein)

        for jobname, probename in self.settings.atomic_probes.items():
            instruction = self.InstructionFile(jobname=jobname,
                                               probename=probename,
                                               settings=self.settings,
                                               cavity=cavity_origin)

            cmds.append('{}'.format(self.settings.superstar_executable)
                        + ' '
                        + '{}.ins'.format(jobname)
                        )

            with PushDir(out):
                instruction.write("{}.ins".format(jobname))
        return cmds

    @staticmethod
    def _merge_cavities(results):
        """
        private method

        atomic hotspot results for the same atomic probe but different cavities are combined onto a single grid
        :param a `hotspots.atomic_hotspot_calculation._AtomicHotspotResult` instance results: a result object
        :return: list, of merged _AtomicHotspotResults
        """
        result_dict = {}
        for result in results:
            if result.identifier in result_dict:
                result_dict[result.identifier].append(result)
            else:
                result_dict.update({result.identifier: [result]})

        merged_results = []
        for identifier, atomic_results in result_dict.items():
            g_dict = {"grid_{}".format(i): g.grid for i, g in enumerate(atomic_results)}
            g = Grid.get_single_grid(g_dict, mask=False)

            b_dict = {"buriedness_{}".format(i): g.buriedness for i, g in enumerate(atomic_results)}
            b = Grid.get_single_grid(b_dict, mask=False)

            merged_results.append(_AtomicHotspotResult(identifier=identifier,
                                                      grid=g,
                                                      buriedness=b,
                                                      )
                                  )

        return merged_results

    def calculate(self, protein, nthreads=None, cavity_origins=None):
        """ Calculates the Atomic Hotspot

        This function executes the Atomic Hotspot Calculation for a given input protein.


        :param `ccdc.protein.Protein` protein: The input protein for which the Atomic Hotspot is to be calculated
        :param int nthreads: The number of processor to be used in the calculation. NB: do not exceed the number of available CPU's
        :param list cavity_origins: The list of cavity origins, if supplied the Atomic Hotspot detection will run on a cavity mode. This increase the speed of the calculation however small cavity can be missed in some cases.

        :return: list of :class:`hotspots._AtomicHotspotResults` instances


        >>> from pdb_python_api import PDBResult
        >>> from ccdc.protein import Protein
        >>> from hotspots.atomic_hotspot_calculation import _AtomicHotspot

        >>> if __name__ == "__main__":
        >>>     # NB: main guard required for multiprocessing on windows!
        >>>     PDBResult("1mtx").download(out_dir="./1mtx")
        >>>     protein = Protein.from_file("1mtx.pdb")
        >>>     protein.add_hydrogens()
        >>>     protein.remove_all_waters()

        >>>     a = _AtomicHotspot()
        >>>     a.settings.atomic_probes = {"apolar": "AROMATIC CH CARBON",
        >>>                                 "donor": "UNCHARGED NH NITROGEN",
        >>>                                 "acceptor": "CARBONYL OXYGEN"}
        >>>     results = a.calculate(protein=protein, nthreads=3)
        [_AtomicHotspotResult(donor), _AtomicHotspotResult(apolar), _AtomicHotspotResult(acceptor)]

        """
        self._merge = False
        if cavity_origins:
            if len(cavity_origins) > 1:
                self._merge = True
            temp_dirs = []
            cmds = []
            env_str = [self.settings.superstar_env] * len(self.settings.atomic_probes) * len(cavity_origins)
            jobnames = list(self.settings.atomic_probes.keys()) * len(cavity_origins)

            for i, cavity_origin in enumerate(cavity_origins):
                out = (join(self.settings.temp_dir, str(i)))
                if not exists(out):
                    mkdir(out)
                t = []
                t.append(out)
                t *= len(self.settings.atomic_probes)
                temp_dirs.extend(t)
                cmds.extend(self._get_cmd(protein, cavity_origin, out=out))

        else:
            temp_dirs = [self.settings.temp_dir] * len(self.settings.atomic_probes)
            env_str = [self.settings.superstar_env] * len(self.settings.atomic_probes)
            jobnames = self.settings.atomic_probes.keys()
            cmds = self._get_cmd(protein, cavity_origins)

        inputs = zip(cmds, jobnames, env_str, temp_dirs)

        results = []
        if nthreads:
            with futures.ProcessPoolExecutor(max_workers=nthreads) as executor:
                for t, j in executor.map(_run_job, inputs):
                    results.append(_AtomicHotspotResult.find(temp_dir=t, jobname=j))

        else:
            for input in inputs:
                t, j = _run_job(input)
                results.append(_AtomicHotspotResult.find(temp_dir=t, jobname=j))

        if self._merge is True:
            results = self._merge_cavities(results)

        return results


class _AtomicHotspotResult(object):
    """
    A class to handle the results of the Atomic Hotspot calculation

    :param str identifier: identifier of the _AtomicHotspotResult, <probe_type>
    :param `hotspots.grid_extension.Grid` grid: atomic propensity grid from Atomic Hotspot calculation
    :param `hotspots.grid_extension.Grid` buriedness: pocket grid from Atomic Hotspot calculation
    :param str ins: the input data to the Atomic Hotspot calculations (superstar ins)
    """

    def __init__(self, identifier, grid, buriedness, ins=None):
        self._identifier = identifier
        self._grid = grid
        self._buriedness = buriedness
        self._ins = ins

    def __str__(self):
        return '_AtomicHotspotResult({})'.format(self.identifier)
    __repr__ = __str__

    @property
    def identifier(self):
        """
        identifier of the _AtomicHotspotResult, <probe_type>
        :return:
        """
        return self._identifier

    @property
    def grid(self):
        """
        atomic propensity grid from Atomic Hotspot calculation
        :return:
        """
        return self._grid

    @grid.setter
    def grid(self, g):
        """
        atomic propensity grid from Atomic Hotspot calculation
        :param `hotspots.grid_extension.Grid` g: a grid
        :return:
        """
        if isinstance(g, Grid):
            self._grid = g
        else:
            raise IOError("Must supply a hotspots.grid_extension.Grid class instance")

    @property
    def buriedness(self):
        """
        pocket grid from Atomic Hotspot calculation
        :return:
        """
        return self._buriedness

    @staticmethod
    def find(temp_dir, jobname):
        """
        searches the calculation working directory and constructs a
        `hotspots.atomic_hotspots_calculation.AtomHotspotResult` instance
        :param str temp_dir: path to the temporary calculation directory
        :param str jobname: name of the atomic probe used in the calculation
        :return: a `hotspots.atomic_hotspots_calculation.AtomHotspotResult` instance
        """
        grid_path = join(temp_dir, jobname + ".acnt")
        grid = Grid.from_file(grid_path)
        buriedness_path = join(temp_dir, jobname + ".ligsite.acnt")

        if exists(buriedness_path):
            b_grid = Grid.from_file(buriedness_path)
            print('ligsite extrema', b_grid.extrema)
            # b.write('buriedness_broken.grd')
            buriedness = _AtomicHotspotResult._correct_ligsite(grid, b_grid)

        else:
            print(buriedness_path)
            raise AttributeError('{} ligsite grid could not be found'.format(jobname))

        return _AtomicHotspotResult(identifier=jobname,
                                   grid=grid,
                                   buriedness=buriedness,
                                   ins=join(temp_dir, "{}.ins".format(jobname))
                                   )

    @staticmethod
    def _correct_ligsite(g, l):
        """
        a correction applied to the pocket grid created by LIGSITE during the Atomic Hotspot calculation
        grid points where ligsite has a score of 0 (i.e. a clash) and SuperStar has a favourable score, set the ligsite
        grid point to its maximum neighbour
        :param `ccdc.utilities.Grid` g: atomic hotspot grid
        :param `ccdc.utilities.Grid` l: pocket grid
        :return:
        """
        mask = ((l < 1) & (g > 2))
        lc = l.copy()
        lc = lc.max_value_of_neighbours().max_value_of_neighbours()
        return l + ((mask * lc).mean_value_of_neighbours() * mask)


def _atomic_hotspot_ins(jobname, probename, settings):
    """
    template for atomic hotspot input file (SuperStar ins)
    :param str jobname: general format (<probe_type>.ins)
    :param str probename: probe identifier
    :param `hotspots.atomic_hotspot_calculation.AtomicHotspot.Settings`settings:
    :return: str, atomic hotspot calculation input str
    """
    ss_str = '''
JOBNAME {0}
PROBENAME {1}
GRID_SPACING 0.5
CALC_CONTOUR_SURFACE 0
SAVE_SUPER_SCATTERPLOT 0
RUN_TYPE SUPERSTAR
DATA_SOURCE CSD
SHELL_VALUES 0.5 0
SIGCHECK_METHOD POISSON SHELL
PROPENSITY_CORRECTION LOGP DEFAULT DEFAULT
MIN_CAVITY_VOLUME 10
CAVITY_RADIUS 10
PEAK_FITTING 0
PEAK_FITTING_NCYCLES 1
MIN_PEAK_HEIGHT 0
PEAK_FITTING_REFINE 0
MOLECULE_FILE protein.pdb
CAVITY_DETECTION 1
MIN_PSP 5
SAVE_CAVITY MESH
USE_TORSIONS 1
RMS_ACCEPT 0.1
NO_MAP_WHEN_INCOMPLETE 0
MAP_FORMAT SYBYL_ASCII_CONTOUR
FIELD_TYPE PROPENSITY
MAP_DISTRIBUTE_POINTS GAUSSIAN_SMEARING
MAX_PENETRATION 0.3
BOX_BORDER {3}
SUPERSTAR_MAP_BACKGROUND {2}
MULTIPLY_POLAR_ONLY 0
INVERSE_MAP 0
AVOID_CLOSE_GRIDPOINTS 0
MAX_NOISE_LEVEL 1
SIGNIFICANCE_LEVEL 0.001
MIN_PROPENSITY {4}
SAVE_SUPERSTAR_MAP 1
SAVE_CENTRAL_GROUPS 0
SAVE_EXCLUDED_VOLUME_MAP 0
SAVE_ACCESSIBLE_VOLUME_MAP 0
FLEXIBILITY_EXHAUSTIVE 0
FLEXIBILITY_FLAG FLAG_ENSEMBLE
SMEARING_SIGMA {5}
SMEARING_NSIGMAS 2
POLAR_POLAR_CORRECTION 1
POLAR_APOLAR_CORRECTION 1
APOLAR_APOLAR_CORRECTION 9
REMOVE_BACKBONE_CONTACTS 1
REMOVE_SUSPICIOUS_CONTACTS 0
MAX_HBOND_DISTANCE 3.1
MIN_HBOND_ANGLE 100
PROPER_HBOND_DISTANCE 2.9
PROPER_HBOND_ANGLE 140
SAVE_LIGSITE_MAP 1
MIN_PEAK_DISTANCE 1
SAVE_PEAKS 1
SAVE_PEAK_MAP 0
CALC_METAL_BONDS 1
SAVE_GOLD_FITTING_POINTS 0
GOLD_MIN_PROPENSITY 5
MAX_SCATTER_ATOMS 5000
MIN_SCATTER_DENS 1.01
SURFACE_FORMAT SYBYL_MOL2
CALC_CONNOLLY_SURFACE 0
CONTOUR_LIST 2 CYAN 4 PURPLE 8 ORANGE
COVALENT_TOLERANCE 0.4
SCORE_PACKING_SHELL 0
ATOM_SCORING SINGLE_POINT'''.format(jobname,
                                    probename,
                                    settings.mapbackgroundvalue,
                                    settings.boxborder,
                                    settings.minpropensity,
                                    settings.superstar_sigma)

    return ss_str
