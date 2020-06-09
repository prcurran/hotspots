
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

_mysphere = [COLOR, 0.9411764705882353, 0.20392156862745098, 0.20392156862745098] +  [ALPHA, 0.5] + [SPHERE, float(0), float(0), float(0), float(2)]

cmd.load_cgo(_mysphere, "mysphere", 1)
_mysphere1 = [COLOR, 0.9411764705882353, 0.20392156862745098, 0.20392156862745098] +  [ALPHA, 0.5] + [SPHERE, float(1), float(1), float(1), float(1)]

cmd.load_cgo(_mysphere1, "mysphere1", 1)