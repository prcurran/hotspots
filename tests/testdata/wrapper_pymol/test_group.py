
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

cmd.pseudoatom(object="myline1pa1", pos=(0, 0, 0), color=(1, 1, 1))

cmd.pseudoatom(object="myline1pa2", pos=(1, 1, 1), color=(1, 1, 1))

cmd.distance(name="myline1", selection1="myline1pa1", selection2="myline1pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.9411764705882353, 0.20392156862745098, 0.20392156862745098), selection="myline1")
cmd.set("dash_width", 10)
cmd.delete("myline1pa1")
cmd.delete("myline1pa2")
cmd.pseudoatom(object="mypseudoatom1", pos=(0, 0, 0), color=(0.9411764705882353, 0.20392156862745098, 0.20392156862745098, 0.5))

cmd.pseudoatom(object="mypseudoatom2", pos=(1, 1, 1), color=(0.9411764705882353, 0.20392156862745098, 0.20392156862745098, 0.5))

cmd.group("mygroup", members="myline1")
cmd.group("mygroup", members="mypseudoatom1")
cmd.group("mygroup", members="mypseudoatom2")