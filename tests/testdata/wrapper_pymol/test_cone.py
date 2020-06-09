
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

_mycone = [CONE, 0, 0, 0, 10, 10, 10, 5, 0, 1, 0, 0, 1, 0, 0, 1, 0.0] 

cmd.load_cgo(_mycone, "mycone", 1)