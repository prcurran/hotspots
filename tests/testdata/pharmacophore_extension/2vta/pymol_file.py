
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

dirpath = None
    
cmd.load("ligands.mol2", "ligands")
cmd.load("protein.mol2", "protein")
obj_donor_projected_point_0 =[COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(27.308), float(3.547), float(65.078), float(0.2)]
cmd.load_cgo(obj_donor_projected_point_0, "donor_projected_point_0", 1)
obj_donor_projected_projection_0 =[COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(26.877), float(0.821), float(65.99), float(0.2)]
cmd.load_cgo(obj_donor_projected_projection_0, "donor_projected_projection_0", 1)
cmd.pseudoatom(object="donor_projected_line_0pa1", pos=[27.308, 3.547, 65.078], color=(1, 1, 1))

cmd.pseudoatom(object="donor_projected_line_0pa2", pos=[26.877, 0.821, 65.99], color=(1, 1, 1))

cmd.distance(name="donor_projected_line_0", selection1="donor_projected_line_0pa1", selection2="donor_projected_line_0pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_projected_line_0")
cmd.set("dash_width", 4.0)
cmd.delete("donor_projected_line_0pa1")
cmd.delete("donor_projected_line_0pa2")
cmd.select("sele", "resi 81")
cmd.show("sticks", "sele")
cmd.group("donor_projected_0", members="donor_projected_point_0")
cmd.group("donor_projected_0", members="donor_projected_projection_0")
cmd.group("donor_projected_0", members="donor_projected_line_0")
cmd.group("donor_projected_pts", members="donor_projected_0")
obj_donor_ch_projected_point_1 =[COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(28.384), float(5.413), float(65.29), float(0.2)]
cmd.load_cgo(obj_donor_ch_projected_point_1, "donor_ch_projected_point_1", 1)
obj_donor_ch_projected_projection_1 =[COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(31.078), float(7.686), float(65.613), float(0.2)]
cmd.load_cgo(obj_donor_ch_projected_projection_1, "donor_ch_projected_projection_1", 1)
cmd.pseudoatom(object="donor_ch_projected_line_1pa1", pos=[28.384, 5.413, 65.29], color=(1, 1, 1))

cmd.pseudoatom(object="donor_ch_projected_line_1pa2", pos=[31.078, 7.686, 65.613], color=(1, 1, 1))

cmd.distance(name="donor_ch_projected_line_1", selection1="donor_ch_projected_line_1pa1", selection2="donor_ch_projected_line_1pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_ch_projected_line_1")
cmd.set("dash_width", 4.0)
cmd.delete("donor_ch_projected_line_1pa1")
cmd.delete("donor_ch_projected_line_1pa2")
cmd.select("sele", "resi 2043")
cmd.show("sticks", "sele")
cmd.group("donor_ch_projected_1", members="donor_ch_projected_point_1")
cmd.group("donor_ch_projected_1", members="donor_ch_projected_projection_1")
cmd.group("donor_ch_projected_1", members="donor_ch_projected_line_1")
cmd.group("donor_ch_projected_pts", members="donor_ch_projected_1")
obj_acceptor_projected_point_2 =[COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(28.364), float(4.137), float(65.698), float(0.2)]
cmd.load_cgo(obj_acceptor_projected_point_2, "acceptor_projected_point_2", 1)
obj_acceptor_projected_projection_2 =[COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(30.481), float(2.301), float(66.766), float(0.2)]
cmd.load_cgo(obj_acceptor_projected_projection_2, "acceptor_projected_projection_2", 1)
cmd.pseudoatom(object="acceptor_projected_line_2pa1", pos=[28.364, 4.137, 65.698], color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_2pa2", pos=[30.481, 2.301, 66.766], color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_2", selection1="acceptor_projected_line_2pa1", selection2="acceptor_projected_line_2pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_2")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_2pa1")
cmd.delete("acceptor_projected_line_2pa2")
cmd.select("sele", "resi 83")
cmd.show("sticks", "sele")
cmd.group("acceptor_projected_2", members="acceptor_projected_point_2")
cmd.group("acceptor_projected_2", members="acceptor_projected_projection_2")
cmd.group("acceptor_projected_2", members="acceptor_projected_line_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
obj_ring_point_3 =[COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(27.606399999999997), float(4.635400000000001), float(64.9408), float(0.2)]
cmd.load_cgo(obj_ring_point_3, "ring_point_3", 1)
obj_ring_projection_3 =[COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(25.021), float(4.551), float(67.359), float(0.2)]
cmd.load_cgo(obj_ring_projection_3, "ring_projection_3", 1)
cmd.pseudoatom(object="ring_line_3pa1", pos=[27.606399999999997, 4.635400000000001, 64.9408], color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_3pa2", pos=[25.021, 4.551, 67.359], color=(1, 1, 1))

cmd.distance(name="ring_line_3", selection1="ring_line_3pa1", selection2="ring_line_3pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_3")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_3pa1")
cmd.delete("ring_line_3pa2")
cmd.select("sele", "resi 31")
cmd.show("sticks", "sele")
cmd.group("ring_3", members="ring_point_3")
cmd.group("ring_3", members="ring_projection_3")
cmd.group("ring_3", members="ring_line_3")
cmd.group("ring_pts", members="ring_3")
obj_ring_point_4 =[COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(27.606399999999997), float(4.635400000000001), float(64.9408), float(0.2)]
cmd.load_cgo(obj_ring_point_4, "ring_point_4", 1)
obj_ring_projection_4 =[COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(29.212), float(3.588), float(62.356), float(0.2)]
cmd.load_cgo(obj_ring_projection_4, "ring_projection_4", 1)
cmd.pseudoatom(object="ring_line_4pa1", pos=[27.606399999999997, 4.635400000000001, 64.9408], color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_4pa2", pos=[29.212, 3.588, 62.356], color=(1, 1, 1))

cmd.distance(name="ring_line_4", selection1="ring_line_4pa1", selection2="ring_line_4pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_4")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_4pa1")
cmd.delete("ring_line_4pa2")
cmd.select("sele", "resi 134")
cmd.show("sticks", "sele")
cmd.group("ring_4", members="ring_point_4")
cmd.group("ring_4", members="ring_projection_4")
cmd.group("ring_4", members="ring_line_4")
cmd.group("ring_pts", members="ring_3")
cmd.group("ring_pts", members="ring_4")
obj_ring_point_5 =[COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(27.606399999999997), float(4.635400000000001), float(64.9408), float(0.2)]
cmd.load_cgo(obj_ring_point_5, "ring_point_5", 1)
obj_ring_projection_5 =[COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(29.97783333333333), float(4.491333333333333), float(69.818), float(0.2)]
cmd.load_cgo(obj_ring_projection_5, "ring_projection_5", 1)
cmd.pseudoatom(object="ring_line_5pa1", pos=[27.606399999999997, 4.635400000000001, 64.9408], color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_5pa2", pos=[29.97783333333333, 4.491333333333333, 69.818], color=(1, 1, 1))

cmd.distance(name="ring_line_5", selection1="ring_line_5pa1", selection2="ring_line_5pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_5")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_5pa1")
cmd.delete("ring_line_5pa2")
cmd.select("sele", "resi 82")
cmd.show("sticks", "sele")
cmd.group("ring_5", members="ring_point_5")
cmd.group("ring_5", members="ring_projection_5")
cmd.group("ring_5", members="ring_line_5")
cmd.group("ring_pts", members="ring_3")
cmd.group("ring_pts", members="ring_4")
cmd.group("ring_pts", members="ring_5")
obj_ring_point_6 =[COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(26.208333333333332), float(5.523333333333333), float(63.53716666666667), float(0.2)]
cmd.load_cgo(obj_ring_point_6, "ring_point_6", 1)
obj_ring_projection_6 =[COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(23.626), float(7.911), float(65.964), float(0.2)]
cmd.load_cgo(obj_ring_projection_6, "ring_projection_6", 1)
cmd.pseudoatom(object="ring_line_6pa1", pos=[26.208333333333332, 5.523333333333333, 63.53716666666667], color=(1, 1, 1))

cmd.pseudoatom(object="ring_line_6pa2", pos=[23.626, 7.911, 65.964], color=(1, 1, 1))

cmd.distance(name="ring_line_6", selection1="ring_line_6pa1", selection2="ring_line_6pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_line_6")
cmd.set("dash_width", 4.0)
cmd.delete("ring_line_6pa1")
cmd.delete("ring_line_6pa2")
cmd.select("sele", "resi 18")
cmd.show("sticks", "sele")
cmd.group("ring_6", members="ring_point_6")
cmd.group("ring_6", members="ring_projection_6")
cmd.group("ring_6", members="ring_line_6")
cmd.group("ring_pts", members="ring_3")
cmd.group("ring_pts", members="ring_4")
cmd.group("ring_pts", members="ring_5")
cmd.group("ring_pts", members="ring_6")
cmd.group("ligand_pharmacophore", members="donor_projected_pts")
cmd.group("ligand_pharmacophore", members="donor_ch_projected_pts")
cmd.group("ligand_pharmacophore", members="acceptor_projected_pts")
cmd.group("ligand_pharmacophore", members="ring_pts")