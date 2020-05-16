
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

dirpath = None
    
cmd.load("ligands.mol2", "ligands")
obj_donor_projected_point_0 =[COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(29.3), float(42.9), float(34.8), float(0.5)]
cmd.load_cgo(obj_donor_projected_point_0, "donor_projected_point_0", 1)
obj_donor_projected_projection_0 =[COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(26.7), float(43.7), float(34.0), float(0.5)]
cmd.load_cgo(obj_donor_projected_projection_0, "donor_projected_projection_0", 1)
cmd.pseudoatom(object="donor_projected_line_0pa1", pos=[29.3, 42.9, 34.8], color=(1, 1, 1))

cmd.pseudoatom(object="donor_projected_line_0pa2", pos=[26.7, 43.7, 34.0], color=(1, 1, 1))

cmd.distance(name="donor_projected_line_0", selection1="donor_projected_line_0pa1", selection2="donor_projected_line_0pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_projected_line_0")
cmd.set("dash_width", 4.0)
cmd.delete("donor_projected_line_0pa1")
cmd.delete("donor_projected_line_0pa2")
cmd.group("donor_projected_0", members="donor_projected_point_0")
cmd.group("donor_projected_0", members="donor_projected_projection_0")
cmd.group("donor_projected_0", members="donor_projected_line_0")
cmd.group("donor_projected_pts", members="donor_projected_0")
obj_donor_projected_point_1 =[COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 1.0] + [SPHERE, float(26.4), float(43.0), float(31.0), float(0.5)]
cmd.load_cgo(obj_donor_projected_point_1, "donor_projected_point_1", 1)
cmd.group("donor_projected_1", members="donor_projected_point_1")
cmd.group("donor_projected_pts", members="donor_projected_0")
cmd.group("donor_projected_pts", members="donor_projected_1")
obj_acceptor_projected_point_2 =[COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(27.6), float(40.8), float(35.0), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_point_2, "acceptor_projected_point_2", 1)
obj_acceptor_projected_projection_2 =[COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(27.3), float(43.2), float(33.5), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_projection_2, "acceptor_projected_projection_2", 1)
cmd.pseudoatom(object="acceptor_projected_line_2pa1", pos=[27.6, 40.8, 35.0], color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_2pa2", pos=[27.3, 43.2, 33.5], color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_2", selection1="acceptor_projected_line_2pa1", selection2="acceptor_projected_line_2pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_2")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_2pa1")
cmd.delete("acceptor_projected_line_2pa2")
cmd.group("acceptor_projected_2", members="acceptor_projected_point_2")
cmd.group("acceptor_projected_2", members="acceptor_projected_projection_2")
cmd.group("acceptor_projected_2", members="acceptor_projected_line_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
obj_acceptor_projected_point_3 =[COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(28.8), float(39.6), float(36.4), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_point_3, "acceptor_projected_point_3", 1)
obj_acceptor_projected_projection_3 =[COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(26.9), float(37.6), float(36.3), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_projection_3, "acceptor_projected_projection_3", 1)
cmd.pseudoatom(object="acceptor_projected_line_3pa1", pos=[28.8, 39.6, 36.4], color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_3pa2", pos=[26.9, 37.6, 36.3], color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_3", selection1="acceptor_projected_line_3pa1", selection2="acceptor_projected_line_3pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_3")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_3pa1")
cmd.delete("acceptor_projected_line_3pa2")
cmd.group("acceptor_projected_3", members="acceptor_projected_point_3")
cmd.group("acceptor_projected_3", members="acceptor_projected_projection_3")
cmd.group("acceptor_projected_3", members="acceptor_projected_line_3")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_3")
obj_acceptor_projected_point_4 =[COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(31.5), float(43.1), float(34.6), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_point_4, "acceptor_projected_point_4", 1)
obj_acceptor_projected_projection_4 =[COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(33.8), float(44.4), float(33.4), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_projection_4, "acceptor_projected_projection_4", 1)
cmd.pseudoatom(object="acceptor_projected_line_4pa1", pos=[31.5, 43.1, 34.6], color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_4pa2", pos=[33.8, 44.4, 33.4], color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_4", selection1="acceptor_projected_line_4pa1", selection2="acceptor_projected_line_4pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_4")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_4pa1")
cmd.delete("acceptor_projected_line_4pa2")
cmd.group("acceptor_projected_4", members="acceptor_projected_point_4")
cmd.group("acceptor_projected_4", members="acceptor_projected_projection_4")
cmd.group("acceptor_projected_4", members="acceptor_projected_line_4")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_3")
cmd.group("acceptor_projected_pts", members="acceptor_projected_4")
obj_acceptor_projected_point_5 =[COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(31.5), float(43.1), float(34.6), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_point_5, "acceptor_projected_point_5", 1)
obj_acceptor_projected_projection_5 =[COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(31.9), float(41.0), float(36.2), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_projection_5, "acceptor_projected_projection_5", 1)
cmd.pseudoatom(object="acceptor_projected_line_5pa1", pos=[31.5, 43.1, 34.6], color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_5pa2", pos=[31.9, 41.0, 36.2], color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_5", selection1="acceptor_projected_line_5pa1", selection2="acceptor_projected_line_5pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_5")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_5pa1")
cmd.delete("acceptor_projected_line_5pa2")
cmd.group("acceptor_projected_5", members="acceptor_projected_point_5")
cmd.group("acceptor_projected_5", members="acceptor_projected_projection_5")
cmd.group("acceptor_projected_5", members="acceptor_projected_line_5")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_3")
cmd.group("acceptor_projected_pts", members="acceptor_projected_4")
cmd.group("acceptor_projected_pts", members="acceptor_projected_5")
obj_acceptor_projected_point_6 =[COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(30.2), float(42.6), float(31.2), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_point_6, "acceptor_projected_point_6", 1)
cmd.group("acceptor_projected_6", members="acceptor_projected_point_6")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_3")
cmd.group("acceptor_projected_pts", members="acceptor_projected_4")
cmd.group("acceptor_projected_pts", members="acceptor_projected_5")
cmd.group("acceptor_projected_pts", members="acceptor_projected_6")
obj_acceptor_projected_point_7 =[COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(29.0), float(44.3), float(29.8), float(0.5)]
cmd.load_cgo(obj_acceptor_projected_point_7, "acceptor_projected_point_7", 1)
cmd.group("acceptor_projected_7", members="acceptor_projected_point_7")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_3")
cmd.group("acceptor_projected_pts", members="acceptor_projected_4")
cmd.group("acceptor_projected_pts", members="acceptor_projected_5")
cmd.group("acceptor_projected_pts", members="acceptor_projected_6")
cmd.group("acceptor_projected_pts", members="acceptor_projected_7")
obj_ring_planar_projected_point_8 =[COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(32.4), float(46.7), float(35.6), float(0.5)]
cmd.load_cgo(obj_ring_planar_projected_point_8, "ring_planar_projected_point_8", 1)
obj_ring_planar_projected_projection_8 =[COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(33.1), float(44.0), float(35.9), float(0.5)]
cmd.load_cgo(obj_ring_planar_projected_projection_8, "ring_planar_projected_projection_8", 1)
cmd.pseudoatom(object="ring_planar_projected_line_8pa1", pos=[32.4, 46.7, 35.6], color=(1, 1, 1))

cmd.pseudoatom(object="ring_planar_projected_line_8pa2", pos=[33.1, 44.0, 35.9], color=(1, 1, 1))

cmd.distance(name="ring_planar_projected_line_8", selection1="ring_planar_projected_line_8pa1", selection2="ring_planar_projected_line_8pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_planar_projected_line_8")
cmd.set("dash_width", 4.0)
cmd.delete("ring_planar_projected_line_8pa1")
cmd.delete("ring_planar_projected_line_8pa2")
cmd.group("ring_planar_projected_8", members="ring_planar_projected_point_8")
cmd.group("ring_planar_projected_8", members="ring_planar_projected_projection_8")
cmd.group("ring_planar_projected_8", members="ring_planar_projected_line_8")
cmd.group("ring_planar_projected_pts", members="ring_planar_projected_8")
obj_ring_planar_projected_point_9 =[COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 1.0] + [SPHERE, float(32.4), float(46.7), float(35.6), float(0.5)]
cmd.load_cgo(obj_ring_planar_projected_point_9, "ring_planar_projected_point_9", 1)
obj_ring_planar_projected_projection_9 =[COLOR, 0.73, 1.0, 0.6] +  [ALPHA, 1.0] + [SPHERE, float(31.7), float(49.4), float(35.2), float(0.5)]
cmd.load_cgo(obj_ring_planar_projected_projection_9, "ring_planar_projected_projection_9", 1)
cmd.pseudoatom(object="ring_planar_projected_line_9pa1", pos=[32.4, 46.7, 35.6], color=(1, 1, 1))

cmd.pseudoatom(object="ring_planar_projected_line_9pa2", pos=[31.7, 49.4, 35.2], color=(1, 1, 1))

cmd.distance(name="ring_planar_projected_line_9", selection1="ring_planar_projected_line_9pa1", selection2="ring_planar_projected_line_9pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.73, 1.0, 0.6), selection="ring_planar_projected_line_9")
cmd.set("dash_width", 4.0)
cmd.delete("ring_planar_projected_line_9pa1")
cmd.delete("ring_planar_projected_line_9pa2")
cmd.group("ring_planar_projected_9", members="ring_planar_projected_point_9")
cmd.group("ring_planar_projected_9", members="ring_planar_projected_projection_9")
cmd.group("ring_planar_projected_9", members="ring_planar_projected_line_9")
cmd.group("ring_planar_projected_pts", members="ring_planar_projected_8")
cmd.group("ring_planar_projected_pts", members="ring_planar_projected_9")
cmd.group("ligand_pharmacophore", members="donor_projected_pts")
cmd.group("ligand_pharmacophore", members="acceptor_projected_pts")
cmd.group("ligand_pharmacophore", members="ring_planar_projected_pts")
cmd.set_color("donor_projected_color", (0.0, 0.0, 1.0))
cmd.load("donor_projected.grd", "donor_projected", state=1)
cmd.isosurface(name="surface_donor_projected", map="donor_projected", level="1")

cmd.color("donor_projected_color", "surface_donor_projected")
cmd.set_color("acceptor_projected_color", (1.0, 0.0, 0.0))
cmd.load("acceptor_projected.grd", "acceptor_projected", state=1)
cmd.isosurface(name="surface_acceptor_projected", map="acceptor_projected", level="1")

cmd.color("acceptor_projected_color", "surface_acceptor_projected")
cmd.set_color("ring_planar_projected_color", (0.33, 1.0, 0.0))
cmd.load("ring_planar_projected.grd", "ring_planar_projected", state=1)
cmd.isosurface(name="surface_ring_planar_projected", map="ring_planar_projected", level="1")

cmd.color("ring_planar_projected_color", "surface_ring_planar_projected")
cmd.group("grds", members="donor_projected")
cmd.group("grds", members="acceptor_projected")
cmd.group("grds", members="ring_planar_projected")
cmd.group("grds", members="surface_donor_projected")
cmd.group("grds", members="surface_acceptor_projected")
cmd.group("grds", members="surface_ring_planar_projected")