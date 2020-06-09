
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

cmd.pseudoatom(object="mylinepa1", pos=(1, 1, 1), color=(1, 1, 1))

cmd.pseudoatom(object="mylinepa2", pos=(2, 2, 2), color=(1, 1, 1))

cmd.distance(name="myline", selection1="mylinepa1", selection2="mylinepa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.9411764705882353, 0.20392156862745098, 0.20392156862745098), selection="myline")
cmd.set("dash_width", 4)
cmd.delete("mylinepa1")
cmd.delete("mylinepa2")
cmd.pseudoatom(object="myline2pa1", pos=(1, 1, 1), color=(1, 1, 1))

cmd.pseudoatom(object="myline2pa2", pos=(0, 0, 0), color=(1, 1, 1))

cmd.distance(name="myline2", selection1="myline2pa1", selection2="myline2pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.9411764705882353, 0.20392156862745098, 0.20392156862745098), selection="myline2")
cmd.set("dash_width", 10)
cmd.delete("myline2pa1")
cmd.delete("myline2pa2")