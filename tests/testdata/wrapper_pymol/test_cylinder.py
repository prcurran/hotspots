
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

_mycylinder = [CYLINDER, 0, 0, 0, 1, 1, 1, 2, 1, 0, 0, 1, 0, 0] 

cmd.load_cgo(_mycylinder, "mycylinder", 1)
_mycylinder2 = [CYLINDER, 0, 0, 0, 10, 10, 10, 0.2, 1, 0, 0, 0.5, 0.5, 0] 

cmd.load_cgo(_mycylinder2, "mycylinder2", 1)