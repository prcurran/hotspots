
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
cmd.load("hotspot_1/apolar.grd", "apolar_hotspot_1")
cmd.isosurface(name="surface_apolar_hotspot_1", map="apolar_hotspot_1", level="5")

cmd.color("yellow", "surface_apolar_hotspot_1")
cmd.set("transparency", 0.2, "surface_apolar_hotspot_1")
cmd.load("hotspot_1/donor.grd", "donor_hotspot_1")
cmd.isosurface(name="surface_donor_hotspot_1", map="donor_hotspot_1", level="5")

cmd.color("blue", "surface_donor_hotspot_1")
cmd.set("transparency", 0.2, "surface_donor_hotspot_1")
cmd.load("hotspot_1/acceptor.grd", "acceptor_hotspot_1")
cmd.isosurface(name="surface_acceptor_hotspot_1", map="acceptor_hotspot_1", level="5")

cmd.color("red", "surface_acceptor_hotspot_1")
cmd.set("transparency", 0.2, "surface_acceptor_hotspot_1")
cmd.group("hotspot_1", members="apolar_hotspot_1")
cmd.group("hotspot_1", members="donor_hotspot_1")
cmd.group("hotspot_1", members="acceptor_hotspot_1")
cmd.group("hotspot_1", members="surface_apolar_hotspot_1")
cmd.group("hotspot_1", members="surface_donor_hotspot_1")
cmd.group("hotspot_1", members="surface_acceptor_hotspot_1")
cmd.pseudoatom(object="PS_apolar_hotspot_1_0", pos=(3.0, 16.5, -29.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_1", pos=(2.5, 17.5, -30.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_2", pos=(1.6666666666666714, 13.25, -34.0), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_3", pos=(2.0, 16.0, -30.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_4", pos=(0.75, 18.333333333333336, -18.75), color=(1, 1, 1), label=23.1)

cmd.pseudoatom(object="PS_apolar_hotspot_1_5", pos=(0.25, 16.0, -35.25), color=(1, 1, 1), label=15.2)

cmd.pseudoatom(object="PS_apolar_hotspot_1_6", pos=(0.5, 10.0, -36.0), color=(1, 1, 1), label=18.9)

cmd.pseudoatom(object="PS_apolar_hotspot_1_7", pos=(-1.0, 14.0, -31.333333333333336), color=(1, 1, 1), label=17.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_8", pos=(-1.0, 14.0, -21.0), color=(1, 1, 1), label=17.6)

cmd.pseudoatom(object="PS_apolar_hotspot_1_9", pos=(-1.5, 9.0, -37.5), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_apolar_hotspot_1_10", pos=(-1.5, 10.0, -37.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_11", pos=(-2.0, 14.0, -34.5), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_apolar_hotspot_1_12", pos=(-3.5, 15.0, -25.5), color=(1, 1, 1), label=30.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_13", pos=(-4.0, 16.0, -24.5), color=(1, 1, 1), label=30.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_14", pos=(-4.5, 16.0, -26.0), color=(1, 1, 1), label=30.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_15", pos=(-17.5, -7.5, -26.0), color=(1, 1, 1), label=14.5)

cmd.pseudoatom(object="PS_apolar_hotspot_1_16", pos=(-23.0, -2.0, -31.0), color=(1, 1, 1), label=14.3)

cmd.pseudoatom(object="PS_apolar_hotspot_1_17", pos=(-22.5, -6.5, -31.5), color=(1, 1, 1), label=23.9)

cmd.pseudoatom(object="PS_apolar_hotspot_1_18", pos=(-23.5, -7.5, -31.0), color=(1, 1, 1), label=20.4)

cmd.pseudoatom(object="PS_apolar_hotspot_1_19", pos=(-23.5, -3.0, -22.0), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_apolar_hotspot_1_20", pos=(-23.5, -3.0, -30.0), color=(1, 1, 1), label=18.6)

cmd.pseudoatom(object="PS_apolar_hotspot_1_21", pos=(-24.25, -7.75, -22.0), color=(1, 1, 1), label=19.9)

cmd.pseudoatom(object="PS_apolar_hotspot_1_22", pos=(-25.0, -3.0, -24.25), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_apolar_hotspot_1_23", pos=(-25.5, -2.0, -24.5), color=(1, 1, 1), label=15.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_24", pos=(-28.0, -3.0, -27.5), color=(1, 1, 1), label=15.5)

cmd.pseudoatom(object="PS_apolar_hotspot_1_25", pos=(-28.5, -4.0, -26.5), color=(1, 1, 1), label=11.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_26", pos=(-29.0, -7.0, -24.5), color=(1, 1, 1), label=19.1)

cmd.pseudoatom(object="PS_apolar_hotspot_1_27", pos=(-29.0, -7.5, -26.0), color=(1, 1, 1), label=19.1)

cmd.pseudoatom(object="PS_apolar_hotspot_1_28", pos=(-30.0, -8.0, -25.0), color=(1, 1, 1), label=19.1)

cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_0")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_0")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_1")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_1")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_2")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_2")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_3")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_3")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_4")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_4")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_5")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_5")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_6")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_6")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_7")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_7")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_8")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_8")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_9")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_9")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_10")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_10")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_11")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_11")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_12")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_12")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_13")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_13")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_14")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_14")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_15")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_15")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_16")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_16")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_17")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_17")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_18")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_18")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_19")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_19")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_20")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_20")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_21")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_21")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_22")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_22")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_23")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_23")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_24")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_24")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_25")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_25")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_26")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_26")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_27")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_27")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_28")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_28")
cmd.pseudoatom(object="PS_donor_hotspot_1_0", pos=(4.5, 19.5, -18.5), color=(1, 1, 1), label=10.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_1", pos=(2.5, 15.0, -34.0), color=(1, 1, 1), label=14.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_2", pos=(1.5, 15.0, -30.0), color=(1, 1, 1), label=12.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_3", pos=(0.5, 9.0, -37.0), color=(1, 1, 1), label=10.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_4", pos=(0.0, 19.0, -19.0), color=(1, 1, 1), label=17.4)

cmd.pseudoatom(object="PS_donor_hotspot_1_5", pos=(-0.5, 16.5, -29.0), color=(1, 1, 1), label=9.6)

cmd.pseudoatom(object="PS_donor_hotspot_1_6", pos=(-0.5, 15.0, -26.0), color=(1, 1, 1), label=9.1)

cmd.pseudoatom(object="PS_donor_hotspot_1_7", pos=(-0.5, 13.5, -30.0), color=(1, 1, 1), label=9.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_8", pos=(-0.5, 12.5, -36.5), color=(1, 1, 1), label=6.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_9", pos=(-1.0, 16.5, -19.5), color=(1, 1, 1), label=20.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_10", pos=(-2.0, 14.0, -22.0), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_11", pos=(-4.0, 16.0, -24.0), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_12", pos=(-19.5, -5.0, -28.5), color=(1, 1, 1), label=7.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_13", pos=(-20.5, -5.5, -34.5), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_14", pos=(-21.0, -5.0, -23.5), color=(1, 1, 1), label=8.9)

cmd.pseudoatom(object="PS_donor_hotspot_1_15", pos=(-22.0, -2.0, -32.5), color=(1, 1, 1), label=10.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_16", pos=(-23.5, -3.0, -31.5), color=(1, 1, 1), label=11.5)

cmd.pseudoatom(object="PS_donor_hotspot_1_17", pos=(-23.5, -8.0, -22.5), color=(1, 1, 1), label=10.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_18", pos=(-23.5, -8.0, -31.0), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_donor_hotspot_1_19", pos=(-26.0, -3.0, -24.5), color=(1, 1, 1), label=10.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_20", pos=(-28.0, -7.0, -27.5), color=(1, 1, 1), label=18.1)

cmd.pseudoatom(object="PS_donor_hotspot_1_21", pos=(-29.5, -6.5, -24.0), color=(1, 1, 1), label=18.6)

cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_0")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_0")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_1")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_1")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_2")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_2")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_3")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_3")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_4")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_4")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_5")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_5")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_6")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_6")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_7")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_7")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_8")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_8")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_9")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_9")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_10")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_10")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_11")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_11")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_12")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_12")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_13")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_13")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_14")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_14")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_15")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_15")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_16")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_16")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_17")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_17")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_18")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_18")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_19")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_19")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_20")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_20")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_21")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_21")
cmd.pseudoatom(object="PS_acceptor_hotspot_1_0", pos=(4.5, 19.5, -18.5), color=(1, 1, 1), label=13.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_1", pos=(2.5, 15.0, -29.0), color=(1, 1, 1), label=19.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_2", pos=(2.5, 15.0, -34.0), color=(1, 1, 1), label=19.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_3", pos=(0.5, 15.0, -15.5), color=(1, 1, 1), label=13.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_4", pos=(0.5, 9.0, -37.0), color=(1, 1, 1), label=13.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_5", pos=(0.0, 19.0, -19.0), color=(1, 1, 1), label=25.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_6", pos=(0.0, 12.5, -36.0), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_7", pos=(-0.5, 16.5, -19.0), color=(1, 1, 1), label=34.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_8", pos=(-0.5, 14.5, -26.0), color=(1, 1, 1), label=11.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_9", pos=(-0.5, 13.5, -29.5), color=(1, 1, 1), label=14.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_10", pos=(-1.0, 16.5, -29.0), color=(1, 1, 1), label=14.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_11", pos=(-1.5, 13.5, -22.0), color=(1, 1, 1), label=18.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_12", pos=(-2.0, 18.5, -16.5), color=(1, 1, 1), label=21.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_13", pos=(-4.0, 16.0, -24.0), color=(1, 1, 1), label=14.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_14", pos=(-16.5, -4.5, -25.0), color=(1, 1, 1), label=12.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_15", pos=(-19.5, -5.0, -28.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_16", pos=(-20.5, -5.5, -34.5), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_17", pos=(-21.0, -5.0, -23.5), color=(1, 1, 1), label=17.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_18", pos=(-22.0, -1.0, -31.5), color=(1, 1, 1), label=14.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_19", pos=(-23.5, -3.0, -31.5), color=(1, 1, 1), label=15.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_20", pos=(-23.5, -8.0, -22.5), color=(1, 1, 1), label=15.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_21", pos=(-23.5, -8.0, -31.0), color=(1, 1, 1), label=14.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_22", pos=(-26.0, -2.5, -24.5), color=(1, 1, 1), label=13.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_23", pos=(-26.5, -5.0, -21.0), color=(1, 1, 1), label=12.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_24", pos=(-28.0, -7.0, -27.5), color=(1, 1, 1), label=20.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_25", pos=(-29.5, -6.5, -24.0), color=(1, 1, 1), label=18.4)

cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_0")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_0")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_1")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_1")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_2")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_2")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_3")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_3")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_4")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_4")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_5")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_5")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_6")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_6")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_7")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_7")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_8")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_8")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_9")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_9")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_10")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_10")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_11")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_11")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_12")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_12")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_13")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_13")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_14")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_14")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_15")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_15")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_16")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_16")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_17")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_17")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_18")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_18")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_19")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_19")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_20")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_20")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_21")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_21")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_22")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_22")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_23")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_23")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_24")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_24")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_25")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_25")
cmd.group("labels_hotspot_1", members="label_apolar_hotspot_1")
cmd.group("labels_hotspot_1", members="label_donor_hotspot_1")
cmd.group("labels_hotspot_1", members="label_acceptor_hotspot_1")
cmd.load("hotspot_1/protein.pdb", "protein_hotspot_1")
cmd.show("cartoon", "protein_hotspot_1")
cmd.hide("line", "protein_hotspot_1")
cmd.show("sticks", "organic")
cmd.load("hotspot_1/apolar.grd", "apolar_hotspot_1")
cmd.isosurface(name="surface_apolar_hotspot_1", map="apolar_hotspot_1", level="5")

cmd.color("yellow", "surface_apolar_hotspot_1")
cmd.set("transparency", 0.2, "surface_apolar_hotspot_1")
cmd.load("hotspot_1/donor.grd", "donor_hotspot_1")
cmd.isosurface(name="surface_donor_hotspot_1", map="donor_hotspot_1", level="5")

cmd.color("blue", "surface_donor_hotspot_1")
cmd.set("transparency", 0.2, "surface_donor_hotspot_1")
cmd.load("hotspot_1/acceptor.grd", "acceptor_hotspot_1")
cmd.isosurface(name="surface_acceptor_hotspot_1", map="acceptor_hotspot_1", level="5")

cmd.color("red", "surface_acceptor_hotspot_1")
cmd.set("transparency", 0.2, "surface_acceptor_hotspot_1")
cmd.group("hotspot_1", members="apolar_hotspot_1")
cmd.group("hotspot_1", members="donor_hotspot_1")
cmd.group("hotspot_1", members="acceptor_hotspot_1")
cmd.group("hotspot_1", members="surface_apolar_hotspot_1")
cmd.group("hotspot_1", members="surface_donor_hotspot_1")
cmd.group("hotspot_1", members="surface_acceptor_hotspot_1")
cmd.pseudoatom(object="PS_apolar_hotspot_1_0", pos=(3.0, 16.5, -29.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_1", pos=(2.5, 17.5, -30.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_2", pos=(1.6666666666666714, 13.25, -34.0), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_3", pos=(2.0, 16.0, -30.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_4", pos=(0.75, 18.333333333333336, -18.75), color=(1, 1, 1), label=23.1)

cmd.pseudoatom(object="PS_apolar_hotspot_1_5", pos=(0.25, 16.0, -35.25), color=(1, 1, 1), label=15.2)

cmd.pseudoatom(object="PS_apolar_hotspot_1_6", pos=(0.5, 10.0, -36.0), color=(1, 1, 1), label=18.9)

cmd.pseudoatom(object="PS_apolar_hotspot_1_7", pos=(-1.0, 14.0, -31.333333333333336), color=(1, 1, 1), label=17.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_8", pos=(-1.0, 14.0, -21.0), color=(1, 1, 1), label=17.6)

cmd.pseudoatom(object="PS_apolar_hotspot_1_9", pos=(-1.5, 9.0, -37.5), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_apolar_hotspot_1_10", pos=(-1.5, 10.0, -37.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_1_11", pos=(-2.0, 14.0, -34.5), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_apolar_hotspot_1_12", pos=(-3.5, 15.0, -25.5), color=(1, 1, 1), label=30.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_13", pos=(-4.0, 16.0, -24.5), color=(1, 1, 1), label=30.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_14", pos=(-4.5, 16.0, -26.0), color=(1, 1, 1), label=30.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_15", pos=(-17.5, -7.5, -26.0), color=(1, 1, 1), label=14.5)

cmd.pseudoatom(object="PS_apolar_hotspot_1_16", pos=(-23.0, -2.0, -31.0), color=(1, 1, 1), label=14.3)

cmd.pseudoatom(object="PS_apolar_hotspot_1_17", pos=(-22.5, -6.5, -31.5), color=(1, 1, 1), label=23.9)

cmd.pseudoatom(object="PS_apolar_hotspot_1_18", pos=(-23.5, -7.5, -31.0), color=(1, 1, 1), label=20.4)

cmd.pseudoatom(object="PS_apolar_hotspot_1_19", pos=(-23.5, -3.0, -22.0), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_apolar_hotspot_1_20", pos=(-23.5, -3.0, -30.0), color=(1, 1, 1), label=18.6)

cmd.pseudoatom(object="PS_apolar_hotspot_1_21", pos=(-24.25, -7.75, -22.0), color=(1, 1, 1), label=19.9)

cmd.pseudoatom(object="PS_apolar_hotspot_1_22", pos=(-25.0, -3.0, -24.25), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_apolar_hotspot_1_23", pos=(-25.5, -2.0, -24.5), color=(1, 1, 1), label=15.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_24", pos=(-28.0, -3.0, -27.5), color=(1, 1, 1), label=15.5)

cmd.pseudoatom(object="PS_apolar_hotspot_1_25", pos=(-28.5, -4.0, -26.5), color=(1, 1, 1), label=11.8)

cmd.pseudoatom(object="PS_apolar_hotspot_1_26", pos=(-29.0, -7.0, -24.5), color=(1, 1, 1), label=19.1)

cmd.pseudoatom(object="PS_apolar_hotspot_1_27", pos=(-29.0, -7.5, -26.0), color=(1, 1, 1), label=19.1)

cmd.pseudoatom(object="PS_apolar_hotspot_1_28", pos=(-30.0, -8.0, -25.0), color=(1, 1, 1), label=19.1)

cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_0")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_0")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_1")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_1")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_2")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_2")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_3")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_3")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_4")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_4")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_5")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_5")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_6")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_6")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_7")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_7")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_8")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_8")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_9")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_9")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_10")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_10")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_11")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_11")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_12")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_12")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_13")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_13")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_14")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_14")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_15")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_15")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_16")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_16")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_17")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_17")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_18")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_18")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_19")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_19")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_20")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_20")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_21")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_21")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_22")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_22")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_23")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_23")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_24")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_24")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_25")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_25")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_26")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_26")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_27")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_27")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_28")
cmd.group("label_apolar_hotspot_1", members="PS_apolar_hotspot_1_28")
cmd.pseudoatom(object="PS_donor_hotspot_1_0", pos=(4.5, 19.5, -18.5), color=(1, 1, 1), label=10.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_1", pos=(2.5, 15.0, -34.0), color=(1, 1, 1), label=14.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_2", pos=(1.5, 15.0, -30.0), color=(1, 1, 1), label=12.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_3", pos=(0.5, 9.0, -37.0), color=(1, 1, 1), label=10.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_4", pos=(0.0, 19.0, -19.0), color=(1, 1, 1), label=17.4)

cmd.pseudoatom(object="PS_donor_hotspot_1_5", pos=(-0.5, 16.5, -29.0), color=(1, 1, 1), label=9.6)

cmd.pseudoatom(object="PS_donor_hotspot_1_6", pos=(-0.5, 15.0, -26.0), color=(1, 1, 1), label=9.1)

cmd.pseudoatom(object="PS_donor_hotspot_1_7", pos=(-0.5, 13.5, -30.0), color=(1, 1, 1), label=9.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_8", pos=(-0.5, 12.5, -36.5), color=(1, 1, 1), label=6.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_9", pos=(-1.0, 16.5, -19.5), color=(1, 1, 1), label=20.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_10", pos=(-2.0, 14.0, -22.0), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_11", pos=(-4.0, 16.0, -24.0), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_12", pos=(-19.5, -5.0, -28.5), color=(1, 1, 1), label=7.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_13", pos=(-20.5, -5.5, -34.5), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_14", pos=(-21.0, -5.0, -23.5), color=(1, 1, 1), label=8.9)

cmd.pseudoatom(object="PS_donor_hotspot_1_15", pos=(-22.0, -2.0, -32.5), color=(1, 1, 1), label=10.2)

cmd.pseudoatom(object="PS_donor_hotspot_1_16", pos=(-23.5, -3.0, -31.5), color=(1, 1, 1), label=11.5)

cmd.pseudoatom(object="PS_donor_hotspot_1_17", pos=(-23.5, -8.0, -22.5), color=(1, 1, 1), label=10.3)

cmd.pseudoatom(object="PS_donor_hotspot_1_18", pos=(-23.5, -8.0, -31.0), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_donor_hotspot_1_19", pos=(-26.0, -3.0, -24.5), color=(1, 1, 1), label=10.0)

cmd.pseudoatom(object="PS_donor_hotspot_1_20", pos=(-28.0, -7.0, -27.5), color=(1, 1, 1), label=18.1)

cmd.pseudoatom(object="PS_donor_hotspot_1_21", pos=(-29.5, -6.5, -24.0), color=(1, 1, 1), label=18.6)

cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_0")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_0")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_1")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_1")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_2")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_2")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_3")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_3")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_4")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_4")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_5")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_5")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_6")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_6")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_7")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_7")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_8")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_8")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_9")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_9")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_10")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_10")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_11")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_11")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_12")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_12")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_13")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_13")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_14")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_14")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_15")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_15")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_16")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_16")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_17")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_17")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_18")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_18")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_19")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_19")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_20")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_20")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_21")
cmd.group("label_donor_hotspot_1", members="PS_donor_hotspot_1_21")
cmd.pseudoatom(object="PS_acceptor_hotspot_1_0", pos=(4.5, 19.5, -18.5), color=(1, 1, 1), label=13.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_1", pos=(2.5, 15.0, -29.0), color=(1, 1, 1), label=19.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_2", pos=(2.5, 15.0, -34.0), color=(1, 1, 1), label=19.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_3", pos=(0.5, 15.0, -15.5), color=(1, 1, 1), label=13.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_4", pos=(0.5, 9.0, -37.0), color=(1, 1, 1), label=13.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_5", pos=(0.0, 19.0, -19.0), color=(1, 1, 1), label=25.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_6", pos=(0.0, 12.5, -36.0), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_7", pos=(-0.5, 16.5, -19.0), color=(1, 1, 1), label=34.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_8", pos=(-0.5, 14.5, -26.0), color=(1, 1, 1), label=11.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_9", pos=(-0.5, 13.5, -29.5), color=(1, 1, 1), label=14.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_10", pos=(-1.0, 16.5, -29.0), color=(1, 1, 1), label=14.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_11", pos=(-1.5, 13.5, -22.0), color=(1, 1, 1), label=18.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_12", pos=(-2.0, 18.5, -16.5), color=(1, 1, 1), label=21.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_13", pos=(-4.0, 16.0, -24.0), color=(1, 1, 1), label=14.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_14", pos=(-16.5, -4.5, -25.0), color=(1, 1, 1), label=12.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_15", pos=(-19.5, -5.0, -28.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_16", pos=(-20.5, -5.5, -34.5), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_17", pos=(-21.0, -5.0, -23.5), color=(1, 1, 1), label=17.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_18", pos=(-22.0, -1.0, -31.5), color=(1, 1, 1), label=14.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_19", pos=(-23.5, -3.0, -31.5), color=(1, 1, 1), label=15.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_20", pos=(-23.5, -8.0, -22.5), color=(1, 1, 1), label=15.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_21", pos=(-23.5, -8.0, -31.0), color=(1, 1, 1), label=14.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_22", pos=(-26.0, -2.5, -24.5), color=(1, 1, 1), label=13.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_23", pos=(-26.5, -5.0, -21.0), color=(1, 1, 1), label=12.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_24", pos=(-28.0, -7.0, -27.5), color=(1, 1, 1), label=20.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_1_25", pos=(-29.5, -6.5, -24.0), color=(1, 1, 1), label=18.4)

cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_0")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_0")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_1")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_1")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_2")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_2")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_3")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_3")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_4")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_4")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_5")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_5")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_6")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_6")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_7")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_7")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_8")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_8")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_9")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_9")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_10")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_10")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_11")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_11")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_12")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_12")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_13")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_13")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_14")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_14")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_15")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_15")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_16")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_16")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_17")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_17")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_18")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_18")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_19")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_19")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_20")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_20")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_21")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_21")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_22")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_22")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_23")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_23")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_24")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_24")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_25")
cmd.group("label_acceptor_hotspot_1", members="PS_acceptor_hotspot_1_25")
cmd.group("labels_hotspot_1", members="label_apolar_hotspot_1")
cmd.group("labels_hotspot_1", members="label_donor_hotspot_1")
cmd.group("labels_hotspot_1", members="label_acceptor_hotspot_1")
cmd.load("hotspot_1/protein.pdb", "protein_hotspot_1")
cmd.show("cartoon", "protein_hotspot_1")
cmd.hide("line", "protein_hotspot_1")
cmd.show("sticks", "organic")


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


surface_list = {'hotspot_1': ['surface_apolar_hotspot_1', 'surface_donor_hotspot_1', 'surface_acceptor_hotspot_1']}

top = tk.Toplevel(plugins.get_tk_root())
master = tk.Frame(top, padx=10, pady=10)
master.pack(fill="both", expand=1)

for child in list(master.children.values()):
    child.destroy()

mmf = tk.Frame(master)

tk.Label(mmf, text=f'0').grid(row=0, column=1, sticky='w')
tk.Label(mmf, text=f'34.0').grid(row=0, column=3, sticky='e')
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
                     to=34.0, resolution=0.1, showvalue=0, variable=v)
        e.grid(row=2 + k, column=2, sticky="ew")
        e = tk.Entry(master, textvariable=v, width=4)
        e.grid(row=2 + k, column=3, sticky="e")

        master.columnconfigure(2, weight=1)

    j += 1



cmd.bg_color("white")
if wd:
    os.chdir(wd)