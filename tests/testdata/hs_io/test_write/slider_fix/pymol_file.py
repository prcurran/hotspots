
from os.path import join
import tempfile
import tkinter as tk
import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

finish_launching()

dirpath = tempfile.mkdtemp()
zip_dir = "out.zip"
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

os.chdir(dirpath)
cmd.load("hotspot/acceptor.grd", "acceptor_hotspot")
cmd.isosurface(name="surface_acceptor_hotspot", map="acceptor_hotspot", level="5")

cmd.color("red", "surface_acceptor_hotspot")
cmd.set("transparency", 0.2, "surface_acceptor_hotspot")
cmd.load("hotspot/donor.grd", "donor_hotspot")
cmd.isosurface(name="surface_donor_hotspot", map="donor_hotspot", level="5")

cmd.color("blue", "surface_donor_hotspot")
cmd.set("transparency", 0.2, "surface_donor_hotspot")
cmd.load("hotspot/apolar.grd", "apolar_hotspot")
cmd.isosurface(name="surface_apolar_hotspot", map="apolar_hotspot", level="5")

cmd.color("yellow", "surface_apolar_hotspot")
cmd.set("transparency", 0.2, "surface_apolar_hotspot")
cmd.group("hotspot", members="acceptor_hotspot")
cmd.group("hotspot", members="donor_hotspot")
cmd.group("hotspot", members="apolar_hotspot")
cmd.group("hotspot", members="surface_acceptor_hotspot")
cmd.group("hotspot", members="surface_donor_hotspot")
cmd.group("hotspot", members="surface_apolar_hotspot")


class IsoLevel(tk.Variable):
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


surface_list = {'hotspot': ['surface_acceptor_hotspot', 'surface_donor_hotspot', 'surface_apolar_hotspot']}

top = tk.Toplevel(plugins.get_tk_root())
master = tk.Frame(top, padx=10, pady=10)
master.pack(fill="both", expand=1)

for child in list(master.children.values()):
    child.destroy()

mmf = tk.Frame(master)

tk.Label(mmf, text=f'0').grid(row=0, column=1, sticky='w')
tk.Label(mmf, text=f'24.4').grid(row=0, column=3, sticky='e')
mmf.grid(row=0, column=2, sticky='ew')
mmf.columnconfigure(1, weight=1)

j = 0
for identifier, surface_list in surface_list.items():
    for i, surface in enumerate(surface_list):
        probe = surface.split("_")[1]
        v = IsoLevel(master, surface, 5)

        k = i + (j * (len(surface_list) + 1))

        if i == 0:
            tk.Label(master, text=identifier).grid(row=1 + k, column=0, sticky="w")

        tk.Label(master, text=probe).grid(row=2 + k, column=1, sticky="w")
        e = tk.Scale(master, orient=tk.HORIZONTAL, from_=0, 
                     to=24.4, resolution=0.1, showvalue=0, variable=v)
        e.grid(row=2 + k, column=2, sticky="ew")
        e = tk.Entry(master, textvariable=v, width=4)
        e.grid(row=2 + k, column=3, sticky="e")

        master.columnconfigure(2, weight=1)

    j += 1



