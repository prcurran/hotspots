import tempfile
import time
import sys
import glob
import subprocess
import collections

from os.path import join, exists, isfile, dirname
from os import environ, mkdir
from concurrent import futures

from ccdc.utilities import PushDir
from ccdc.io import csd_directory, MoleculeWriter

from pprint import pprint

try:
    from grid_extension import Grid
except ImportError:
    from ccdc.utilities import Grid

Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])


def _run_job(args):
    """
    runs atomic hotspot(superstar) job (paralyzable)
    :param args:
    :return:
    """
    cmd, jobname, superstar_env, temp_dir = args
    env = environ.copy()
    env.update(superstar_env)
    with PushDir(temp_dir):
        subprocess.call(cmd, shell=sys.platform != 'win32', env=env)

    return temp_dir, jobname


class AtomicHotspot(object):
    """
    class to handle SuperStar run
    """

    class Settings(object):
        """
        setting for Superstar run
        """

        def __init__(self):
            self.database = 'CSD'
            self.mapbackgroundvalue = 1
            self.boxborder = 10
            self.minpropensity = 1
            self.superstar_sigma = 0.5
            self.superstar_executable, self.superstar_env = self._set_environment_variables()
            self.temp_dir = tempfile.mkdtemp()

            # TODO: decide which probes should be included
            # TODO: this could be read in from superstar defaults (talk to Richard)
            self._csd_atomic_probes = {}

            # TODO: add PDB probes
            self._pdb_atomic_probes = {}

        @staticmethod
        def _set_environment_variables():
            """
            From the internal API (Thanks Richard!), sets up superstar environment variables
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
                    merc = glob.glob(join(base, 'mercury*'))
                    if len(merc):
                        merc = merc[0]

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
        """
        A class to handle 'SuperStar' instruction files
        """

        def __init__(self, jobname, probename, settings, cavity=None):
            self.ins_str = atomic_hotspot_ins(jobname, probename, settings)
            if cavity:
                self.ins_str += '\nCAVITY_ORIGIN {} {} {}'.format(cavity[0], cavity[1], cavity[2])
            else:
                self.ins_str += '\nSUBSTRUCTURE ALL'

        def write(self, fname):
            with open(fname, "w") as f:
                f.write(self.ins_str)

    def __init__(self, settings=None):
        if not settings:
            self.settings = self.Settings()

        else:
            settings = settings

    # @staticmethod
    # def _run_job(args):
    #     """
    #     runs atomic hotspot(superstar) job (paralyzable)
    #     :param args:
    #     :return:
    #     """
    #     cmd, jobname, superstar_env, temp_dir = args
    #     env = environ.copy()
    #     env.update(superstar_env)
    #     with PushDir(temp_dir):
    #         subprocess.call(cmd, shell=sys.platform != 'win32', env=env)
    #     return AtomicHotspotResult.find(temp_dir=temp_dir, jobname=jobname)

    def _get_cmd(self, protein, cavity_origin, out=None):
        """
        constructs the commandline str required to run atomic hotspot calculation
        :param jobname:
        :param probename:
        :param out:
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

            cmds.append('{}'.format(self.settings.superstar_executable) + ' ' + '{}.ins'.format(jobname))
            with PushDir(out):
                instruction.write("{}.ins".format(jobname))
        return cmds

    def _merge_cavities(self, results):
        """
        atomic hotspot results for the same atomic probe but different cavities are combined onto a single grid
        :param results:
        :return:
        """
        # TODO: supply ins list during AtomicHotspotResult creation, (only if requested)
        # TODO: deal with background value of 1 during additions
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


            merged_results.append(AtomicHotspotResult(identifier=identifier,
                                                      grid=g,
                                                      buriedness=b,
                                                      ins = None
                                                      )
                                  )

        return merged_results

    def calculate(self, protein, nthreads=None, cavity_origins=None):
        """
        main function of AtomicHotspot, used to handle atomic hotspot calculations
        :param protein:
        :param probes: MUST BE LIST
        :param nthreads:
        :param cavity_origins: MUST BE LIST
        :return:
        """
        self._merge = False
        if cavity_origins:

            if len(cavity_origins) > 1:
                self._merge = True

            # create input lists per cavity
            temp_dirs = []
            cmds = []
            env_str = [self.settings.superstar_env] * len(self.settings.atomic_probes) * len(cavity_origins)
            jobnames = self.settings.atomic_probes.keys() * len(cavity_origins)

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
            # create input lists
            temp_dirs = [self.settings.temp_dir] * len(self.settings.atomic_probes)
            env_str = [self.settings.superstar_env] * len(self.settings.atomic_probes)
            jobnames = self.settings.atomic_probes.keys()
            cmds = self._get_cmd(protein, cavity_origins)

        inputs = zip(cmds, jobnames, env_str, temp_dirs)

        # paralyze atomic hotspot calculation
        results = []
        if nthreads:
            with futures.ProcessPoolExecutor(max_workers=nthreads) as executor:
                for t, j  in executor.map(_run_job, inputs):
                    results.append(AtomicHotspotResult.find(temp_dir=t, jobname=j) )

        else:
            for input in inputs:
                _run_job(input)
            results = [AtomicHotspotResult.find(temp_dir=t, jobname=j) for t, j in result_input]

        # merge the atomic hotspot results for different cavities
        if self._merge == True:
            results = self._merge_cavities(results)

        return results


class AtomicHotspotResult(object):
    """
    class to store a SuperStar result
    """

    def __init__(self, identifier, grid, buriedness, ins=None):
        """
        attributes of a AtomicHotspotResult object
        :param identifier:
        :param grid:
        :param buriedness:
        """
        self._identifier = identifier
        self._grid = grid
        self._buriedness = buriedness
        self._ins = ins

    @property
    def identifier(self):
        return self._identifier

    @property
    def buriedness(self):
        return self._buriedness

    @property
    def grid(self):
        return self._grid

    @grid.setter
    def grid(self, g):
        if isinstance(g, Grid):
            self._grid = g
        else:
            raise IOError("Must supply a hotspots.grid_extension.Grid class instance")

    @staticmethod
    def find(temp_dir, jobname):
        """
        constructs as AtomicHotspotResult object from AtomicHotspot calculation
        :param settings:
        :return:
        """
        grid_path = join(temp_dir, jobname + ".acnt")

        grid = Grid.from_file(grid_path)

        # if exists(grid_path):
        #     grid = Grid.from_file(grid_path)
        #
        # else:
        #     print grid_path
        #     raise AttributeError('{} superstar grid could not be found'.format(jobname))

        buriedness_path = join(temp_dir, jobname + ".ligsite.acnt")
        if exists(buriedness_path):
            b = Grid.from_file(buriedness_path)
            buriedness = AtomicHotspotResult._correct_ligsite(grid, b)

        else:
            print buriedness_path
            raise AttributeError('{} ligsite grid could not be found'.format(jobname))

        return AtomicHotspotResult(identifier=jobname,
                                   grid=grid,
                                   buriedness=buriedness,
                                   ins = join(temp_dir, "{}.ins".format(jobname)))

    @staticmethod
    def _correct_ligsite(g, l):
        """
        Grid points where ligsite has a score of 0 (i.e. a clash) and SuperStar has a favourable score, set the ligsite
        grid point to its maximum neighbour
        :param g:
        :param l:
        :return:
        """
        mask = ((l < 1) & (g > 2))
        lc = l.copy()
        lc = lc.max_value_of_neighbours().max_value_of_neighbours()
        return l + ((mask * lc).mean_value_of_neighbours() * mask)

def atomic_hotspot_ins(jobname, probename, settings):
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

