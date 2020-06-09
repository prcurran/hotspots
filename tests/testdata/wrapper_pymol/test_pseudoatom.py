
from os.path import join
import tempfile
import tkinter as tk
import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

finish_launching()

cmd.pseudoatom(object="mypseudoatom", pos=(1, 1, 1), color=(0.9411764705882353, 0.20392156862745098, 0.20392156862745098, 0.5), label=None)

cmd.pseudoatom(object="mypseudoatom2", pos=(2, 2, 2), color=(0.9411764705882353, 0.20392156862745098, 0.20392156862745098, 0.5), label=31.2)
