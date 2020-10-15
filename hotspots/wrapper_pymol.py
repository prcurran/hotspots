"""
This module is designed to interface with the PyMOL API.

The PyMOL python package can be installed within conda however we have found launch
new sessions unreliable. Therefore, this module aids the automated creation of a
python script which can be executed manually within PyMOl GUI. We have found this to
be a more reliable way to visualise with PyMOL

For more information on the PyMOL commands available:
    - https://pymol.org/dokuwiki/doku.php?id=api:cmd:alpha

"""
import numpy as np


class PyMOLFile():
    """
    Holds the basics
    """

    def __init__(self):
        self.commands = PyMOLCommands.imports()

    def write(self, path):
        with open(path, 'w') as w:
            w.write(self.commands)


class PyMOLCommands():
    """
    A class containing the implemented PyMOL commands
    """

    @staticmethod
    def imports():
        """
        Provides the required python package imports

        :return: Python code for visualisation in PyMOL
        :rtype: str
        """
        out_str = """
try:
    import tkinter as tk      
except ImportError:
    import Tkinter as tk
from os.path import join
import tempfile

import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

finish_launching()
"""
        return out_str

    @staticmethod
    def load_cgo(obj, objname, state=1):
        """
        Load CGO object.

        :param obj: PyMOL object float list representation
        :param objname: name of object
        :param state: state in which to load

        :type obj: str
        :type objname: str
        :type state: int
        :return:
        """
        return f'\ncmd.load_cgo({obj}, "{objname}", {state})'

    @staticmethod
    def sphere(objname, rgba, coords, radius=1):
        """
        Creates a cgo sphere.

        :param objname: name of object
        :param rgba: RGBA color, values between 0 - 1.
        :param coords: coordinate of sphere centroid
        :param radius: raduis of sphere

        :type objname: str
        :type rgba: (float, float, float, float)
        :type coords: (float, float, float)
        :type radius: float

        :return: PyMOL command
        :rtype: str
        """
        return \
            f'\n{objname} = ' \
            f'[COLOR, {rgba[0]}, {rgba[1]}, {rgba[2]}] + ' \
            f' [ALPHA, {rgba[3]}] +' \
            f' [SPHERE, float({coords[0]}), float({coords[1]}), float({coords[2]}), float({radius})]\n'

    @staticmethod
    def cylinder(objname, coord1, coord2, radius, rgba, rgba2=None):
        """
        Create a CGO cylinder.

        TODO: ALPHA command does not work
        :param objname: name of object
        :param coord1: first coordinate
        :param coord2: second coordinate
        :param radius: cylinder radius
        :param rgba: first color
        :param rgba2: second color

        :type objname: str
        :type coord1: (float, float, float)
        :type coord2: (float, float, float)
        :type radius: float
        :type rgba: (float, float, float, float)
        :type rgba2: (float, float, float)

        :return: PyMOL command
        :rtype: str
        """
        if not rgba2:
            rgba2 = rgba
        return \
            f'\n{objname} = ' \
            f'[CYLINDER, {coord1[0]}, {coord1[1]}, {coord1[2]}, ' \
            f'{coord2[0]}, {coord2[1]}, {coord2[2]}, {radius}, ' \
            f'{rgba[0]}, {rgba[1]}, {rgba[2]}, ' \
            f'{rgba2[0]}, {rgba2[1]}, {rgba2[2]}] \n'

    @staticmethod
    def cone(objname, coord1, coord2, base_radius, rgba, rgba2=None, tip_radius=0, base_alpha=1):
        if not rgba2:
            rgba2 = rgba
        return \
            f'\n{objname} = ' \
            f'[CONE, {coord1[0]}, {coord1[1]}, {coord1[2]}, ' \
            f'{coord2[0]}, {coord2[1]}, {coord2[2]}, ' \
            f'{base_radius}, {tip_radius}, ' \
            f'{rgba[0]}, {rgba[1]}, {rgba[2]}, ' \
            f'{rgba2[0]}, {rgba2[1]}, {rgba2[2]}, ' \
            f'{base_alpha}, 0.0] \n'

    @staticmethod
    def arrow(objname, coord1, coord2, rgba, rgba2=None, radius=0.1, head_length=0.5, arrow_head_ratio=1.75):
        """
        Create a CGO arrow.

        :param objname: name of object
        :param coord1: coordinates of cylinder end
        :param coord2: coordinates of arrow point
        :param rgba: color 1
        :param rgba2: color 2
        :param radius: radius of cylinder
        :param head_length: length of arrow head
        :param arrow_head_ratio: ratio of cylinder radius to cone-base radius

        :type objname: str
        :type coord1: (float, float, float)
        :type coord2: (float, float, float)
        :type rgba: (float, float, float)
        :type rgba2: (float, float, float)
        :type radius: float
        :type head_length: float
        :type arrow_head_ratio: float

        :return: PyMOL command
        :rtype: str
        """
        if not rgba2:
            rgba2 = rgba

        v = np.array(coord2) - np.array(coord1)
        u = np.divide(v, np.linalg.norm(v))
        coord3 = coord2 - (head_length * u)

        return f'{PyMOLCommands.cylinder(f"{objname}_cylinder", coord1, coord3, radius, rgba, rgba2)}' \
               f'{PyMOLCommands.cone(f"{objname}_cone", coord3, coord2, radius * arrow_head_ratio, rgba2)}' \
               f'\n{objname} = {objname}_cylinder + {objname}_cone \n'

    @staticmethod
    def pseudoatom(objname, coords, color=(1, 1, 1), label=None):
        """
        Creates a pseudoatom.

        TODO: Currently, color assignment is not working for pseudoatoms.

        :param objname: name of object
        :param coords: coordinates of pseudoatom
        :param color: rgba of color

        :type objname: str
        :type coords: (float, float, float)
        :type color: (float, float, float)

        :return: PyMOL command
        :rtype: str
        """
        if not label:
            return f'\ncmd.pseudoatom(object="{objname}", pos={coords}, color={color})\n'
        else:
            return f'\ncmd.pseudoatom(object="{objname}", pos={coords}, color={color}, label={label})\n'

    @staticmethod
    def line(objname, coord1, coord2, rgb=(1, 1, 1), width=4.0):
        """
        Creates a line.

        TODO: Current, line width is not working.

        :param objname: name of object
        :param coord1: coordinate of pseudoatom 1
        :param coord2: coordinate of pseudoatom 2
        :param rgb: rgb color
        :param width: width of line

        :type objname: str
        :type coord1: (float, float, float)
        :type coord2: (float, float, float)
        :type rgb: (float, float, float)
        :type width: float

        :return: PyMOL command
        :rtype: str
        """
        pymol_out = PyMOLCommands.pseudoatom(f"{objname}pa1", coord1)
        pymol_out += PyMOLCommands.pseudoatom(f"{objname}pa2", coord2)
        pymol_out += \
            f'\ncmd.distance(name="{objname}", selection1="{objname}pa1", selection2="{objname}pa2", width=0.5, gap=0.2, label=0, state=1)\n' \
            f'\ncmd.set("dash_color", {(rgb[0], rgb[1], rgb[2])}, selection="{objname}")' \
            f'\ncmd.set("dash_width", {width})' \
            f'\ncmd.delete("{objname}pa1")' \
            f'\ncmd.delete("{objname}pa2")'
        return pymol_out

    @staticmethod
    def isosurface(grd_name, isosurface_name, level, color):
        """
        Loads a grid creates a contoured surface.

        :param grd_name: name of grid object
        :param isosurface_name: name of isosurface object
        :param level: value theshold at which to contour
        :param color: named color, NB: must be a registered color

        :type grd_name: str
        :type isosurface_name: str
        :type level: float
        :type color: str

        :return: PyMOL command
        :rtype: str
        """
        # pymol_out = PyMOLCommands.load(fname, grd_name)
        pymol_out = f'\ncmd.isosurface(name="{isosurface_name}", map="{grd_name}", level="{level}")\n'
        pymol_out += f'\ncmd.color("{color}", "{isosurface_name}")'
        return pymol_out

    @staticmethod
    def load(fname, objname=None):
        """
        Loads any PyMOL object

        :param fname: path to file
        :param objname: name of object

        :type fname: str
        :type objname: str

        :return: PyMOL command
        :rtype: str
        """
        if not objname:
            objname = fname.split(".")[0]
        return f'\ncmd.load("{fname}", "{objname}")'

    @staticmethod
    def group(group_name, members):
        """
        Creates a group of PyMOL objects

        :param group_name: name of group
        :param members: list of object names

        :type group_name: str
        :type members: list

        :return: PyMOL command
        :rtype: str
        """
        pymol_out = ""
        for member in members:
            pymol_out += f'\ncmd.group("{group_name}", members="{member}")'
        return pymol_out

    @staticmethod
    def select(objname, selection):
        return \
    f'\ncmd.select("{objname}", "{selection}")'

    @staticmethod
    def show(representation, selection):
        return \
    f'\ncmd.show("{representation}", "{selection}")'

    @staticmethod
    def hide(representation, selection):
        return \
    f'\ncmd.hide("{representation}", "{selection}")'

    @staticmethod
    def set_color(objname, rgb):
        """
        Create a new registered color.

        :param objname: name of new color
        :param rgb: rgb of color (rgba tolerated)

        :type objname: str
        :type rgb: (float, float, float)

        :return: PyMOL command
        :rtype: str
        """
        return f'\ncmd.set_color("{objname}", {(rgb[0], rgb[1], rgb[2])})'

    @staticmethod
    def color(pymol_color, objname):
        """
        For recolor proteins

        :param pymol_color: a registered pymol color (NOT RGBA)
        :param objname: name of obj to color

        :type pymol_color: str
        :type objname: str

        :return: PyMOL command
        :rtype: str
        """
        return f'\ncmd.color("{pymol_color}", "{objname}")'

    @staticmethod
    def pymol_set(setting_name, value, selection):
        return f'\ncmd.set("{setting_name}", {value}, "{selection}")'

    @staticmethod
    def unzip_dir(container):
        return f'\ndirpath = tempfile.mkdtemp()' \
               f'\nzip_dir = "{container}.zip"' \
               f'\nwd = os.getcwd()' \
               f'\nwith zipfile.ZipFile(zip_dir) as hs_zip:' \
               f'\n    hs_zip.extractall(dirpath)' \
               f'\n\nos.chdir(dirpath)'

    @staticmethod
    def push_to_wd():
        return f'\nif wd:\n    os.chdir(wd)'

    @staticmethod
    def background_color(color):
        return f'cmd.bg_color("{color}")'

    @staticmethod
    def isoslider(surface_dic, surface_value_dic, min_value=0):
        """
        Create a GUI to control the isolevel of isosurface loaded

        Adapted from:
        Isocontour slider ("density slider") plugin for PyMOL
        (c) 2013 Thomas Holder
        License: BSD-2-Clause

        :param surface_dic: a dictionary of hotspot identifier by surface labels
        :param min_value: set to 0
        :param max_value: set to the maximum value across all the grids

        :type surface_dic: dict
        :type min_value: float
        :type max_value: float

        :return: PyMOL command
        :rtype: str
        """
        return \
f"""
\n\nclass IsoLevel(tk.Variable):
    def __init__(self, master, name, level):
        tk.Variable.__init__(self, master, value=level)
        self.name = name
        self.trace('w', self.callback)

    def callback(self, *args):
        cmd.isolevel(self.name, self.get())

    def increment(self, event=None, delta=0.1):
        self.set(round(float(self.get()) + delta, 2))

    def decrement(self, event=None):
        self.increment(None, -0.1)


surface_list = {surface_dic}
surface_max_list = {surface_value_dic}

top = tk.Toplevel(plugins.get_tk_root())

master = tk.Frame(top, padx=10, pady=10)
master.pack(fill="both", expand=1)

for child in list(master.children.values()):
    child.destroy()


row_counter = 0
for identifier, component_dic in surface_list.items():
    # add calculation identifier
    tk.Label(master, text=identifier).grid(row=row_counter, column=0, sticky="w")
    row_counter += 1
    
    for component_id, surfaces in component_dic.items():
        # add collection label, e.g. superstar or hotspot etc.
        tk.Label(master, text=component_id).grid(row=row_counter, column=1, sticky='w')
        row_counter += 1
        
        for i, surface in enumerate(surfaces):
            # add grid type label
            probe = surface.split("_")[-2]
            tk.Label(master, text=probe).grid(row=row_counter, column=2, sticky="w")
            
            # slider code 
            v = IsoLevel(master, surface, 5)
            e = tk.Scale(master, orient=tk.HORIZONTAL, from_={min_value}, to=surface_max_list[identifier][component_id],
                         resolution=0.1, showvalue=0, variable=v)
            e.grid(row=row_counter, column=3, sticky="ew")

            e = tk.Entry(master, textvariable=v, width=4)
            e.grid(row=row_counter, column=4, sticky="e")
            master.columnconfigure(3, weight=1)
            row_counter += 1
\n\n
"""
