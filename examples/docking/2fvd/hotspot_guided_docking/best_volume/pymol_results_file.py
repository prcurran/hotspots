
from os.path import join
import tempfile
import zipfile
from pymol import cmd
from pymol.cgo import *

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

    
dirpath = tempfile.mkdtemp()
zip_dir = 'out.zip'
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

cmd.load(join(dirpath,"protein.pdb"), "protein")
cmd.show("cartoon", "protein")

cmd.set("surface_cavity_mode", 1)
cmd.show("surface", "protein")
cmd.set("surface_trim_factor", 13)
cmd.set('transparency', 0.5, "protein")
    
if dirpath:
    f = join(dirpath, "label_threshold_10.mol2")
else:
    f = "label_threshold_10.mol2"

cmd.load(f, 'label_threshold_10')
cmd.hide('everything', 'label_threshold_10')
cmd.label("label_threshold_10", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


if dirpath:
    f = join(dirpath, "label_threshold_14.mol2")
else:
    f = "label_threshold_14.mol2"

cmd.load(f, 'label_threshold_14')
cmd.hide('everything', 'label_threshold_14')
cmd.label("label_threshold_14", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


if dirpath:
    f = join(dirpath, "label_threshold_17.mol2")
else:
    f = "label_threshold_17.mol2"

cmd.load(f, 'label_threshold_17')
cmd.hide('everything', 'label_threshold_17')
cmd.label("label_threshold_17", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10, 14, 17]
gfiles = ['donor.grd', 'acceptor.grd', 'apolar.grd', 'negative.grd']
grids = ['donor', 'acceptor', 'apolar', 'negative']
num = 0
surf_transparency = 0.2

if dirpath:
    gfiles = [join(dirpath, g) for g in gfiles]

for t in threshold_list:
    for i in range(len(grids)):
        try:
            cmd.load(r'%s'%(gfiles[i]), '%s_%s'%(grids[i], str(num)))
            cmd.isosurface('surface_%s_%s_%s'%(grids[i], t, num), '%s_%s'%(grids[i], num), t)
            cmd.set('transparency', surf_transparency, 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.color(colour_dict['%s'%(grids[i])], 'surface_%s_%s_%s'%(grids[i], t, num))
            cmd.group('threshold_%s'%(t), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s' % (t), members='label_threshold_%s' % (t))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s' % (t))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cluster_dict = {"rank_1":[], "rank_1_arrows":[]}

cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-6.0), float(30.5), float(11.5), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([-6.0,30.5,11.5], [-5.47385567252,29.1543346537,14.0014439437], color="blue red", name="Arrows_rank_1_1")

cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-3.0), float(25.5), float(15.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([-3.0,25.5,15.0], [-4.5058811948,27.6713684973,15.2734842763], color="blue red", name="Arrows_rank_1_2")

cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(-0.5), float(27.5), float(13.5), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([-0.5,27.5,13.5], [-0.337991088292,27.6513689537,15.5554932179], color="blue red", name="Arrows_rank_1_3")

cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(26.5), float(11.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([2.5,26.5,11.0], [2.87992406592,27.3303762792,13.7974374753], color="blue red", name="Arrows_rank_1_4")

cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(2.5), float(29.0), float(6.5), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([2.5,29.0,6.5], [4.88487120209,30.9912927316,5.62317829408], color="blue red", name="Arrows_rank_1_5")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(-2.5), float(32.5), float(11.5), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([-2.5,32.5,11.5], [-2.52093353131,31.9912699106,14.3794559293], color="red blue", name="Arrows_rank_1_6")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(1.5), float(31.0), float(6.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([1.5,31.0,6.0], [3.68290289402,33.3592386914,5.43017217444], color="red blue", name="Arrows_rank_1_7")

cluster_dict["rank_1"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-0.840067027133), float(28.8276913067), float(9.71085325574), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(24.5), float(9.0), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 0.6, 0.1, 0.6] + [ALPHA, 0.6] + [SPHERE, float(-5.0), float(27.0), float(10.5), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 0.6, 0.1, 0.6] + [ALPHA, 0.6] + [SPHERE, float(-1.0), float(24.5), float(13.5), float(1.0)]


cmd.load_cgo(cluster_dict["rank_1"], "Features_rank_1", 1)
cmd.load_cgo(cluster_dict["rank_1_arrows"], "Arrows_rank_1")
cmd.set("transparency", 0.2,"Features_rank_1")
cmd.group("Pharmacophore_rank_1", members="Features_rank_1")
cmd.group("Pharmacophore_rank_1", members="Arrows_rank_1")

if dirpath:
    f = join(dirpath, "label_threshold_rank_1.mol2")
else:
    f = "label_threshold_rank_1.mol2"

cmd.load(f, 'label_threshold_rank_1')
cmd.hide('everything', 'label_threshold_rank_1')
cmd.label("label_threshold_rank_1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_rank_1', members= 'label_threshold_rank_1')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
