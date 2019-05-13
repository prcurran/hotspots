
from os.path import join
import tempfile
import zipfile
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

dirpath = None
    
def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.07, gap=0.0, hlength=-1, hradius=-1, color='blue red', name=''):
    from chempy import cpv
    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 +           [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 +           [1.0, 0.0]
    return obj

    
cluster_dict = {"CDK2":[], "CDK2_arrows":[]}

cluster_dict["CDK2"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(22.0), float(12.5), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(6.0), float(22.5), float(8.5), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(31.0), float(7.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-1.5), float(23.5), float(9.5), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-1.5), float(25.0), float(13.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(3.0), float(28.5), float(6.5), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(27.5), float(12.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(5.0), float(22.0), float(2.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(4.5), float(26.0), float(7.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(0.5), float(30.0), float(7.5), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(25.0), float(11.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(26.5), float(8.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-4.0), float(28.0), float(11.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(22.5), float(14.0), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(2.0), float(30.5), float(6.5), float(1.0)]


cluster_dict["CDK2"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(25.5), float(12.0), float(1.0)]


cmd.load_cgo(cluster_dict["CDK2"], "Features_CDK2", 1)
cmd.load_cgo(cluster_dict["CDK2_arrows"], "Arrows_CDK2")
cmd.set("transparency", 0.2,"Features_CDK2")
cmd.group("Pharmacophore_CDK2", members="Features_CDK2")
cmd.group("Pharmacophore_CDK2", members="Arrows_CDK2")

if dirpath:
    f = join(dirpath, "label_threshold_CDK2.mol2")
else:
    f = "label_threshold_CDK2.mol2"

cmd.load(f, 'label_threshold_CDK2')
cmd.hide('everything', 'label_threshold_CDK2')
cmd.label("label_threshold_CDK2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_CDK2', members= 'label_threshold_CDK2')
