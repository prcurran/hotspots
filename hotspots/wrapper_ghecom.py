import os
import tempfile
from ccdc.utilities import PushDir
from ccdc import io
from hotspots.hs_utilities import Helper


class Ghecom(object):
    """
    Minimal wrapper for pocket detection using Ghecom
    NB: Currently this method is only available to linux users
    Please ensure you have set the following environment variable:
    >>>export GHECOM_EXE='<path_to_ghecom>'
    """
    def __init__(self):
        self.temp = tempfile.mkdtemp()
        self.input = "protein.pdb"
        self.output = "ghecom_out.pdb"

    def run(self, protein, grid_spacing=0.5):
        """
        executes the command line call

        NB: different versions of ghecom have different commandline args

        :param protein: protein
        :param grid_spacing: grid spacing, must be the same used in hotspot calculation

        :type protein: `ccdc.protein.Protein`
        :type grid_spacing: float
        """
        with PushDir(self.temp):
            with io.MoleculeWriter('protein.pdb') as writer:
                writer.write(protein)

                cmd = f"{os.environ['GHECOM_EXE']} -ipdb {self.input} -M M -gw {grid_spacing} -rli 2.5 -rlx 9.5 -opocpdb {self.output}"
                os.system(cmd)

        return os.path.join(self.temp, self.output)

    @staticmethod
    def pdb_to_grid(path, template):
        """
        converts pdb file to grid

        :param path: path to the input PDB
        :param template: empty grid, NB: must have same grid spec as superstar grids
        :type path: str
        :type template: `hotspots.grid_extension.Grid`

        :return: populated grid
        :rtype: `hotspots.grid_extension.Grid`
        """

        lines = Helper.get_lines_from_file(path)
        for line in lines:
            if line.startswith("HETATM"):
                coordinates = (float(line[31:38]), float(line[39:46]), float(line[47:54]))
                rinacc = float(line[61:66])
                i, j, k = template.point_to_indices(coordinates)
                nx, ny, nz = template.nsteps
                if 0 < i < nx and 0 < j < ny and 0 < k < nz:
                    template.set_value(i, j, k, 9.5 - rinacc)
        return template
