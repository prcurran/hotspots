
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

cmd.set_color("pomegranate", (0.9411764705882353, 0.20392156862745098, 0.20392156862745098))
cmd.load("test_grid.grd", "mygrid", state=1)
cmd.isosurface(name="mysurface", map="mygrid", level="4")

cmd.color("pomegranate", "mysurface")