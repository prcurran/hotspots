
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
wd = os.getcwd()
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
cmd.pseudoatom(object="PS_acceptor_0", pos=(-17.0, -4.5, -25.5), color=(1, 1, 1), label=12.7)

cmd.pseudoatom(object="PS_acceptor_1", pos=(-19.5, -5.0, -28.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_acceptor_2", pos=(-20.5, -5.5, -35.0), color=(1, 1, 1), label=16.8)

cmd.pseudoatom(object="PS_acceptor_3", pos=(-21.0, -5.0, -23.5), color=(1, 1, 1), label=18.1)

cmd.pseudoatom(object="PS_acceptor_4", pos=(-22.0, -0.5, -31.5), color=(1, 1, 1), label=15.2)

cmd.pseudoatom(object="PS_acceptor_5", pos=(-23.5, -3.0, -31.5), color=(1, 1, 1), label=15.5)

cmd.pseudoatom(object="PS_acceptor_6", pos=(-23.5, -8.0, -22.5), color=(1, 1, 1), label=15.0)

cmd.pseudoatom(object="PS_acceptor_7", pos=(-23.5, -8.0, -31.0), color=(1, 1, 1), label=14.6)

cmd.pseudoatom(object="PS_acceptor_8", pos=(-26.0, -3.0, -24.5), color=(1, 1, 1), label=13.7)

cmd.pseudoatom(object="PS_acceptor_9", pos=(-26.5, -5.0, -21.0), color=(1, 1, 1), label=12.9)

cmd.pseudoatom(object="PS_acceptor_10", pos=(-28.0, -7.0, -27.5), color=(1, 1, 1), label=20.7)

cmd.pseudoatom(object="PS_acceptor_11", pos=(-29.5, -6.5, -24.0), color=(1, 1, 1), label=18.9)

cmd.group("label_acceptor_hotspot", members="PS_acceptor_0")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_0")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_1")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_1")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_2")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_2")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_3")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_3")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_4")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_4")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_5")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_5")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_6")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_6")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_7")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_7")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_8")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_8")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_9")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_9")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_10")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_10")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_11")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_11")
cmd.pseudoatom(object="PS_donor_0", pos=(-19.5, -5.0, -28.5), color=(1, 1, 1), label=7.1)

cmd.pseudoatom(object="PS_donor_1", pos=(-20.5, -5.5, -34.5), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_donor_2", pos=(-21.0, -5.0, -23.5), color=(1, 1, 1), label=8.9)

cmd.pseudoatom(object="PS_donor_3", pos=(-23.5, -3.0, -31.5), color=(1, 1, 1), label=11.7)

cmd.pseudoatom(object="PS_donor_4", pos=(-23.5, -8.0, -22.5), color=(1, 1, 1), label=10.5)

cmd.pseudoatom(object="PS_donor_5", pos=(-23.5, -8.0, -31.0), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_donor_6", pos=(-26.0, -3.0, -24.5), color=(1, 1, 1), label=10.0)

cmd.pseudoatom(object="PS_donor_7", pos=(-28.0, -7.0, -27.5), color=(1, 1, 1), label=18.1)

cmd.pseudoatom(object="PS_donor_8", pos=(-29.5, -6.5, -24.0), color=(1, 1, 1), label=19.6)

cmd.group("label_donor_hotspot", members="PS_donor_0")
cmd.group("label_donor_hotspot", members="PS_donor_0")
cmd.group("label_donor_hotspot", members="PS_donor_1")
cmd.group("label_donor_hotspot", members="PS_donor_1")
cmd.group("label_donor_hotspot", members="PS_donor_2")
cmd.group("label_donor_hotspot", members="PS_donor_2")
cmd.group("label_donor_hotspot", members="PS_donor_3")
cmd.group("label_donor_hotspot", members="PS_donor_3")
cmd.group("label_donor_hotspot", members="PS_donor_4")
cmd.group("label_donor_hotspot", members="PS_donor_4")
cmd.group("label_donor_hotspot", members="PS_donor_5")
cmd.group("label_donor_hotspot", members="PS_donor_5")
cmd.group("label_donor_hotspot", members="PS_donor_6")
cmd.group("label_donor_hotspot", members="PS_donor_6")
cmd.group("label_donor_hotspot", members="PS_donor_7")
cmd.group("label_donor_hotspot", members="PS_donor_7")
cmd.group("label_donor_hotspot", members="PS_donor_8")
cmd.group("label_donor_hotspot", members="PS_donor_8")
cmd.pseudoatom(object="PS_apolar_0", pos=(-17.5, -7.5, -26.0), color=(1, 1, 1), label=16.5)

cmd.pseudoatom(object="PS_apolar_1", pos=(-23.0, -2.0, -31.0), color=(1, 1, 1), label=18.6)

cmd.pseudoatom(object="PS_apolar_2", pos=(-23.5, -7.0, -31.5), color=(1, 1, 1), label=24.4)

cmd.pseudoatom(object="PS_apolar_3", pos=(-23.5, -3.0, -21.5), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_apolar_4", pos=(-24.0, -8.5, -23.5), color=(1, 1, 1), label=22.6)

cmd.pseudoatom(object="PS_apolar_5", pos=(-25.5, -3.5, -24.5), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_6", pos=(-28.5, -3.5, -27.5), color=(1, 1, 1), label=15.6)

cmd.pseudoatom(object="PS_apolar_7", pos=(-29.5, -7.5, -25.5), color=(1, 1, 1), label=19.1)

cmd.group("label_apolar_hotspot", members="PS_apolar_0")
cmd.group("label_apolar_hotspot", members="PS_apolar_0")
cmd.group("label_apolar_hotspot", members="PS_apolar_1")
cmd.group("label_apolar_hotspot", members="PS_apolar_1")
cmd.group("label_apolar_hotspot", members="PS_apolar_2")
cmd.group("label_apolar_hotspot", members="PS_apolar_2")
cmd.group("label_apolar_hotspot", members="PS_apolar_3")
cmd.group("label_apolar_hotspot", members="PS_apolar_3")
cmd.group("label_apolar_hotspot", members="PS_apolar_4")
cmd.group("label_apolar_hotspot", members="PS_apolar_4")
cmd.group("label_apolar_hotspot", members="PS_apolar_5")
cmd.group("label_apolar_hotspot", members="PS_apolar_5")
cmd.group("label_apolar_hotspot", members="PS_apolar_6")
cmd.group("label_apolar_hotspot", members="PS_apolar_6")
cmd.group("label_apolar_hotspot", members="PS_apolar_7")
cmd.group("label_apolar_hotspot", members="PS_apolar_7")
cmd.group("labels", members="label_acceptor_hotspot")
cmd.group("labels", members="label_donor_hotspot")
cmd.group("labels", members="label_apolar_hotspot")
cmd.load("hotspot/protein.pdb", "protein_hotspot")


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




    cmd.bg_color("white")
    cmd.show("cartoon", "protein")
    cmd.color("slate", "protein")
    cmd.show("sticks", "organic")
    cmd.hide("lines", "protein")
    
if wd:
    os.chdir(wd)