
from os.path import join
import tempfile
import tkinter as tk
import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

finish_launching()

cmd.load("ligands.mol2", "ligands")
cmd.load("protein.mol2", "protein")
donor_ch_projected_point_0 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(22.604), float(34.487), float(31.141), float(0.5)]

cmd.load_cgo(donor_ch_projected_point_0, "donor_ch_projected_point_0_obj", 1)
donor_ch_projected_projection_0 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(20.708), float(32.724), float(29.515), float(0.5)]

cmd.load_cgo(donor_ch_projected_projection_0, "donor_ch_projected_projection_0_obj", 1)
cmd.pseudoatom(object="donor_ch_projected_line_0pa1", pos=(22.604, 34.487, 31.141), color=(1, 1, 1))

cmd.pseudoatom(object="donor_ch_projected_line_0pa2", pos=(20.708, 32.724, 29.515), color=(1, 1, 1))

cmd.distance(name="donor_ch_projected_line_0", selection1="donor_ch_projected_line_0pa1", selection2="donor_ch_projected_line_0pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_ch_projected_line_0")
cmd.set("dash_width", 4.0)
cmd.delete("donor_ch_projected_line_0pa1")
cmd.delete("donor_ch_projected_line_0pa2")
cmd.group("donor_ch_projected_0", members="donor_ch_projected_point_0_obj")
cmd.group("donor_ch_projected_0", members="donor_ch_projected_projection_0")
cmd.group("donor_ch_projected_0", members="donor_ch_projected_line_0")
cmd.select("sele", "resi 800")
cmd.show("sticks", "sele")
donor_ch_projected_point_1 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(22.604), float(34.487), float(31.141), float(0.5)]

cmd.load_cgo(donor_ch_projected_point_1, "donor_ch_projected_point_1_obj", 1)
donor_ch_projected_projection_1 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(20.708), float(32.724), float(29.515), float(0.5)]

cmd.load_cgo(donor_ch_projected_projection_1, "donor_ch_projected_projection_1_obj", 1)
cmd.pseudoatom(object="donor_ch_projected_line_1pa1", pos=(22.604, 34.487, 31.141), color=(1, 1, 1))

cmd.pseudoatom(object="donor_ch_projected_line_1pa2", pos=(20.708, 32.724, 29.515), color=(1, 1, 1))

cmd.distance(name="donor_ch_projected_line_1", selection1="donor_ch_projected_line_1pa1", selection2="donor_ch_projected_line_1pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_ch_projected_line_1")
cmd.set("dash_width", 4.0)
cmd.delete("donor_ch_projected_line_1pa1")
cmd.delete("donor_ch_projected_line_1pa2")
cmd.group("donor_ch_projected_1", members="donor_ch_projected_point_1_obj")
cmd.group("donor_ch_projected_1", members="donor_ch_projected_projection_1")
cmd.group("donor_ch_projected_1", members="donor_ch_projected_line_1")
cmd.select("sele", "resi 800")
cmd.show("sticks", "sele")
donor_ch_projected_point_2 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(22.604), float(34.487), float(31.141), float(0.5)]

cmd.load_cgo(donor_ch_projected_point_2, "donor_ch_projected_point_2_obj", 1)
donor_ch_projected_projection_2 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(20.708), float(32.724), float(29.515), float(0.5)]

cmd.load_cgo(donor_ch_projected_projection_2, "donor_ch_projected_projection_2_obj", 1)
cmd.pseudoatom(object="donor_ch_projected_line_2pa1", pos=(22.604, 34.487, 31.141), color=(1, 1, 1))

cmd.pseudoatom(object="donor_ch_projected_line_2pa2", pos=(20.708, 32.724, 29.515), color=(1, 1, 1))

cmd.distance(name="donor_ch_projected_line_2", selection1="donor_ch_projected_line_2pa1", selection2="donor_ch_projected_line_2pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_ch_projected_line_2")
cmd.set("dash_width", 4.0)
cmd.delete("donor_ch_projected_line_2pa1")
cmd.delete("donor_ch_projected_line_2pa2")
cmd.group("donor_ch_projected_2", members="donor_ch_projected_point_2_obj")
cmd.group("donor_ch_projected_2", members="donor_ch_projected_projection_2")
cmd.group("donor_ch_projected_2", members="donor_ch_projected_line_2")
cmd.select("sele", "resi 800")
cmd.show("sticks", "sele")
donor_ch_projected_point_3 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(12.708), float(33.79), float(38.568), float(0.5)]

cmd.load_cgo(donor_ch_projected_point_3, "donor_ch_projected_point_3_obj", 1)
donor_ch_projected_projection_3 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(9.632), float(34.011), float(39.405), float(0.5)]

cmd.load_cgo(donor_ch_projected_projection_3, "donor_ch_projected_projection_3_obj", 1)
cmd.pseudoatom(object="donor_ch_projected_line_3pa1", pos=(12.708, 33.79, 38.568), color=(1, 1, 1))

cmd.pseudoatom(object="donor_ch_projected_line_3pa2", pos=(9.632, 34.011, 39.405), color=(1, 1, 1))

cmd.distance(name="donor_ch_projected_line_3", selection1="donor_ch_projected_line_3pa1", selection2="donor_ch_projected_line_3pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_ch_projected_line_3")
cmd.set("dash_width", 4.0)
cmd.delete("donor_ch_projected_line_3pa1")
cmd.delete("donor_ch_projected_line_3pa2")
cmd.group("donor_ch_projected_3", members="donor_ch_projected_point_3_obj")
cmd.group("donor_ch_projected_3", members="donor_ch_projected_projection_3")
cmd.group("donor_ch_projected_3", members="donor_ch_projected_line_3")
cmd.select("sele", "resi 791")
cmd.show("sticks", "sele")
donor_ch_projected_point_4 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(22.061), float(31.836), float(32.324), float(0.5)]

cmd.load_cgo(donor_ch_projected_point_4, "donor_ch_projected_point_4_obj", 1)
donor_ch_projected_projection_4 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(20.708), float(32.724), float(29.515), float(0.5)]

cmd.load_cgo(donor_ch_projected_projection_4, "donor_ch_projected_projection_4_obj", 1)
cmd.pseudoatom(object="donor_ch_projected_line_4pa1", pos=(22.061, 31.836, 32.324), color=(1, 1, 1))

cmd.pseudoatom(object="donor_ch_projected_line_4pa2", pos=(20.708, 32.724, 29.515), color=(1, 1, 1))

cmd.distance(name="donor_ch_projected_line_4", selection1="donor_ch_projected_line_4pa1", selection2="donor_ch_projected_line_4pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_ch_projected_line_4")
cmd.set("dash_width", 4.0)
cmd.delete("donor_ch_projected_line_4pa1")
cmd.delete("donor_ch_projected_line_4pa2")
cmd.group("donor_ch_projected_4", members="donor_ch_projected_point_4_obj")
cmd.group("donor_ch_projected_4", members="donor_ch_projected_projection_4")
cmd.group("donor_ch_projected_4", members="donor_ch_projected_line_4")
cmd.select("sele", "resi 800")
cmd.show("sticks", "sele")
donor_ch_projected_point_5 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(22.061), float(31.836), float(32.324), float(0.5)]

cmd.load_cgo(donor_ch_projected_point_5, "donor_ch_projected_point_5_obj", 1)
donor_ch_projected_projection_5 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(20.708), float(32.724), float(29.515), float(0.5)]

cmd.load_cgo(donor_ch_projected_projection_5, "donor_ch_projected_projection_5_obj", 1)
cmd.pseudoatom(object="donor_ch_projected_line_5pa1", pos=(22.061, 31.836, 32.324), color=(1, 1, 1))

cmd.pseudoatom(object="donor_ch_projected_line_5pa2", pos=(20.708, 32.724, 29.515), color=(1, 1, 1))

cmd.distance(name="donor_ch_projected_line_5", selection1="donor_ch_projected_line_5pa1", selection2="donor_ch_projected_line_5pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_ch_projected_line_5")
cmd.set("dash_width", 4.0)
cmd.delete("donor_ch_projected_line_5pa1")
cmd.delete("donor_ch_projected_line_5pa2")
cmd.group("donor_ch_projected_5", members="donor_ch_projected_point_5_obj")
cmd.group("donor_ch_projected_5", members="donor_ch_projected_projection_5")
cmd.group("donor_ch_projected_5", members="donor_ch_projected_line_5")
cmd.select("sele", "resi 800")
cmd.show("sticks", "sele")
acceptor_projected_point_6 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(25.053), float(34.204), float(32.067), float(0.5)]

cmd.load_cgo(acceptor_projected_point_6, "acceptor_projected_point_6_obj", 1)
acceptor_projected_projection_6 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(24.543), float(37.258), float(30.116), float(0.5)]

cmd.load_cgo(acceptor_projected_projection_6, "acceptor_projected_projection_6_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_6pa1", pos=(25.053, 34.204, 32.067), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_6pa2", pos=(24.543, 37.258, 30.116), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_6", selection1="acceptor_projected_line_6pa1", selection2="acceptor_projected_line_6pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_6")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_6pa1")
cmd.delete("acceptor_projected_line_6pa2")
cmd.group("acceptor_projected_6", members="acceptor_projected_point_6_obj")
cmd.group("acceptor_projected_6", members="acceptor_projected_projection_6")
cmd.group("acceptor_projected_6", members="acceptor_projected_line_6")
cmd.select("sele", "resi 799")
cmd.show("sticks", "sele")
acceptor_projected_point_7 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(25.053), float(34.204), float(32.067), float(0.5)]

cmd.load_cgo(acceptor_projected_point_7, "acceptor_projected_point_7_obj", 1)
acceptor_projected_projection_7 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(24.543), float(37.258), float(30.116), float(0.5)]

cmd.load_cgo(acceptor_projected_projection_7, "acceptor_projected_projection_7_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_7pa1", pos=(25.053, 34.204, 32.067), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_7pa2", pos=(24.543, 37.258, 30.116), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_7", selection1="acceptor_projected_line_7pa1", selection2="acceptor_projected_line_7pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_7")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_7pa1")
cmd.delete("acceptor_projected_line_7pa2")
cmd.group("acceptor_projected_7", members="acceptor_projected_point_7_obj")
cmd.group("acceptor_projected_7", members="acceptor_projected_projection_7")
cmd.group("acceptor_projected_7", members="acceptor_projected_line_7")
cmd.select("sele", "resi 799")
cmd.show("sticks", "sele")
acceptor_projected_point_8 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(25.053), float(34.204), float(32.067), float(0.5)]

cmd.load_cgo(acceptor_projected_point_8, "acceptor_projected_point_8_obj", 1)
acceptor_projected_projection_8 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(24.543), float(37.258), float(30.116), float(0.5)]

cmd.load_cgo(acceptor_projected_projection_8, "acceptor_projected_projection_8_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_8pa1", pos=(25.053, 34.204, 32.067), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_8pa2", pos=(24.543, 37.258, 30.116), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_8", selection1="acceptor_projected_line_8pa1", selection2="acceptor_projected_line_8pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_8")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_8pa1")
cmd.delete("acceptor_projected_line_8pa2")
cmd.group("acceptor_projected_8", members="acceptor_projected_point_8_obj")
cmd.group("acceptor_projected_8", members="acceptor_projected_projection_8")
cmd.group("acceptor_projected_8", members="acceptor_projected_line_8")
cmd.select("sele", "resi 799")
cmd.show("sticks", "sele")
acceptor_projected_point_9 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(12.683), float(33.305), float(37.276), float(0.5)]

cmd.load_cgo(acceptor_projected_point_9, "acceptor_projected_point_9_obj", 1)
acceptor_projected_projection_9 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(10.035), float(33.489), float(35.966), float(0.5)]

cmd.load_cgo(acceptor_projected_projection_9, "acceptor_projected_projection_9_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_9pa1", pos=(12.683, 33.305, 37.276), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_9pa2", pos=(10.035, 33.489, 35.966), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_9", selection1="acceptor_projected_line_9pa1", selection2="acceptor_projected_line_9pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_9")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_9pa1")
cmd.delete("acceptor_projected_line_9pa2")
cmd.group("acceptor_projected_9", members="acceptor_projected_point_9_obj")
cmd.group("acceptor_projected_9", members="acceptor_projected_projection_9")
cmd.group("acceptor_projected_9", members="acceptor_projected_line_9")
cmd.select("sele", "resi 793")
cmd.show("sticks", "sele")
ring_point_10 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(18.319799999999997), float(31.344), float(34.1718), float(0.5)]

cmd.load_cgo(ring_point_10, "ring_point_10_obj", 1)
ring_projection_10 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(17.316), float(27.744), float(36.095), float(0.5)]

cmd.load_cgo(ring_projection_10, "ring_projection_10_obj", 1)
cmd.pseudoatom(object="ring_line_10pa1", pos=(18.319799999999997, 31.344, 34.1718), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_10pa2", pos=(17.316, 27.744, 36.095), color=(1, 1, 1))

cmd.distance(name="ring_line_10", selection1="ring_line_10pa1", selection2="ring_line_10pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_10")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_10pa1")
cmd.delete("ring_line_10pa2")
cmd.group("ring_10", members="ring_point_10_obj")
cmd.group("ring_10", members="ring_projection_10")
cmd.group("ring_10", members="ring_line_10")
cmd.select("sele", "resi 718")
cmd.show("sticks", "sele")
ring_point_11 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(15.0105), float(32.47333333333333), float(36.04083333333333), float(0.5)]

cmd.load_cgo(ring_point_11, "ring_point_11_obj", 1)
ring_projection_11 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(15.097), float(28.62), float(36.82), float(0.5)]

cmd.load_cgo(ring_projection_11, "ring_projection_11_obj", 1)
cmd.pseudoatom(object="ring_line_11pa1", pos=(15.0105, 32.47333333333333, 36.04083333333333), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_11pa2", pos=(15.097, 28.62, 36.82), color=(1, 1, 1))

cmd.distance(name="ring_line_11", selection1="ring_line_11pa1", selection2="ring_line_11pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_11")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_11pa1")
cmd.delete("ring_line_11pa2")
cmd.group("ring_11", members="ring_point_11_obj")
cmd.group("ring_11", members="ring_projection_11")
cmd.group("ring_11", members="ring_line_11")
cmd.select("sele", "resi 718")
cmd.show("sticks", "sele")
ring_point_12 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(15.0105), float(32.47333333333333), float(36.04083333333333), float(0.5)]

cmd.load_cgo(ring_point_12, "ring_point_12_obj", 1)
ring_projection_12 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(15.608), float(36.305), float(35.323), float(0.5)]

cmd.load_cgo(ring_projection_12, "ring_projection_12_obj", 1)
cmd.pseudoatom(object="ring_line_12pa1", pos=(15.0105, 32.47333333333333, 36.04083333333333), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_12pa2", pos=(15.608, 36.305, 35.323), color=(1, 1, 1))

cmd.distance(name="ring_line_12", selection1="ring_line_12pa1", selection2="ring_line_12pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_12")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_12pa1")
cmd.delete("ring_line_12pa2")
cmd.group("ring_12", members="ring_point_12_obj")
cmd.group("ring_12", members="ring_projection_12")
cmd.group("ring_12", members="ring_line_12")
cmd.select("sele", "resi 844")
cmd.show("sticks", "sele")
ring_point_13 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(16.385666666666665), float(34.15616666666667), float(42.090833333333336), float(0.5)]

cmd.load_cgo(ring_point_13, "ring_point_13_obj", 1)
ring_projection_13 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(13.647), float(36.166), float(40.814), float(0.5)]

cmd.load_cgo(ring_projection_13, "ring_projection_13_obj", 1)
cmd.pseudoatom(object="ring_line_13pa1", pos=(16.385666666666665, 34.15616666666667, 42.090833333333336), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_13pa2", pos=(13.647, 36.166, 40.814), color=(1, 1, 1))

cmd.distance(name="ring_line_13", selection1="ring_line_13pa1", selection2="ring_line_13pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_13")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_13pa1")
cmd.delete("ring_line_13pa2")
cmd.group("ring_13", members="ring_point_13_obj")
cmd.group("ring_13", members="ring_projection_13")
cmd.group("ring_13", members="ring_line_13")
cmd.select("sele", "resi 4")
cmd.show("sticks", "sele")
ring_point_14 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(16.385666666666665), float(34.15616666666667), float(42.090833333333336), float(0.5)]

cmd.load_cgo(ring_point_14, "ring_point_14_obj", 1)
ring_projection_14 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(19.444), float(32.442), float(43.941), float(0.5)]

cmd.load_cgo(ring_projection_14, "ring_projection_14_obj", 1)
cmd.pseudoatom(object="ring_line_14pa1", pos=(16.385666666666665, 34.15616666666667, 42.090833333333336), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_14pa2", pos=(19.444, 32.442, 43.941), color=(1, 1, 1))

cmd.distance(name="ring_line_14", selection1="ring_line_14pa1", selection2="ring_line_14pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_14")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_14pa1")
cmd.delete("ring_line_14pa2")
cmd.group("ring_14", members="ring_point_14_obj")
cmd.group("ring_14", members="ring_projection_14")
cmd.group("ring_14", members="ring_line_14")
cmd.select("sele", "resi 745")
cmd.show("sticks", "sele")
ring_point_15 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(15.216500000000002), float(38.24733333333333), float(44.83783333333333), float(0.5)]

cmd.load_cgo(ring_point_15, "ring_point_15_obj", 1)
ring_projection_15 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(15.63), float(38.183), float(41.174), float(0.5)]

cmd.load_cgo(ring_projection_15, "ring_projection_15_obj", 1)
cmd.pseudoatom(object="ring_line_15pa1", pos=(15.216500000000002, 38.24733333333333, 44.83783333333333), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_15pa2", pos=(15.63, 38.183, 41.174), color=(1, 1, 1))

cmd.distance(name="ring_line_15", selection1="ring_line_15pa1", selection2="ring_line_15pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_15")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_15pa1")
cmd.delete("ring_line_15pa2")
cmd.group("ring_15", members="ring_point_15_obj")
cmd.group("ring_15", members="ring_projection_15")
cmd.group("ring_15", members="ring_line_15")
cmd.select("sele", "resi 854")
cmd.show("sticks", "sele")
ring_point_16 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(15.216500000000002), float(38.24733333333333), float(44.83783333333333), float(0.5)]

cmd.load_cgo(ring_point_16, "ring_point_16_obj", 1)
ring_projection_16 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(16.933), float(38.522), float(49.479), float(0.5)]

cmd.load_cgo(ring_projection_16, "ring_projection_16_obj", 1)
cmd.pseudoatom(object="ring_line_16pa1", pos=(15.216500000000002, 38.24733333333333, 44.83783333333333), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_16pa2", pos=(16.933, 38.522, 49.479), color=(1, 1, 1))

cmd.distance(name="ring_line_16", selection1="ring_line_16pa1", selection2="ring_line_16pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_16")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_16pa1")
cmd.delete("ring_line_16pa2")
cmd.group("ring_16", members="ring_point_16_obj")
cmd.group("ring_16", members="ring_projection_16")
cmd.group("ring_16", members="ring_line_16")
cmd.select("sele", "resi 766")
cmd.show("sticks", "sele")
ring_point_17 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(15.216500000000002), float(38.24733333333333), float(44.83783333333333), float(0.5)]

cmd.load_cgo(ring_point_17, "ring_point_17_obj", 1)
ring_projection_17 = [COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(15.915833333333333), float(43.11933333333333), float(45.16433333333333), float(0.5)]

cmd.load_cgo(ring_projection_17, "ring_projection_17_obj", 1)
cmd.pseudoatom(object="ring_line_17pa1", pos=(15.216500000000002, 38.24733333333333, 44.83783333333333), color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_17pa2", pos=(15.915833333333333, 43.11933333333333, 45.16433333333333), color=(1, 1, 1))

cmd.distance(name="ring_line_17", selection1="ring_line_17pa1", selection2="ring_line_17pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_17")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_17pa1")
cmd.delete("ring_line_17pa2")
cmd.group("ring_17", members="ring_point_17_obj")
cmd.group("ring_17", members="ring_projection_17")
cmd.group("ring_17", members="ring_line_17")
cmd.select("sele", "resi 856")
cmd.show("sticks", "sele")
cmd.group("donor_ch_projected_pts", members="donor_ch_projected_0")
cmd.group("donor_ch_projected_pts", members="donor_ch_projected_1")
cmd.group("donor_ch_projected_pts", members="donor_ch_projected_2")
cmd.group("donor_ch_projected_pts", members="donor_ch_projected_3")
cmd.group("donor_ch_projected_pts", members="donor_ch_projected_4")
cmd.group("donor_ch_projected_pts", members="donor_ch_projected_5")
cmd.group("acceptor_projected_pts", members="acceptor_projected_6")
cmd.group("acceptor_projected_pts", members="acceptor_projected_7")
cmd.group("acceptor_projected_pts", members="acceptor_projected_8")
cmd.group("acceptor_projected_pts", members="acceptor_projected_9")
cmd.group("ring_pts", members="ring_10")
cmd.group("ring_pts", members="ring_11")
cmd.group("ring_pts", members="ring_12")
cmd.group("ring_pts", members="ring_13")
cmd.group("ring_pts", members="ring_14")
cmd.group("ring_pts", members="ring_15")
cmd.group("ring_pts", members="ring_16")
cmd.group("ring_pts", members="ring_17")
cmd.group("ligand_pharmacophore", members="donor_projected_pts")
cmd.group("ligand_pharmacophore", members="donor_ch_projected_pts")
cmd.group("ligand_pharmacophore", members="acceptor_projected_pts")
cmd.group("ligand_pharmacophore", members="ring_pts")