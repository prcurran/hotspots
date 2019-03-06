
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
zip_dir = 'hotspot_boundaries.zip'
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

cmd.load(join(dirpath,"protein.pdb"), "protein")
cmd.show("cartoon", "protein")

cmd.set("surface_cavity_mode", 1)
cmd.show("surface", "protein")
cmd.set("surface_trim_factor", 13)
cmd.set('transparency', 0.5, "protein")
    
if dirpath:
    f = join(dirpath, "0\label_threshold_17.8.mol2")
else:
    f = "0\label_threshold_17.8.mol2"

cmd.load(f, 'label_threshold_17.8')
cmd.hide('everything', 'label_threshold_17.8')
cmd.label("label_threshold_17.8", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [17.8]
gfiles = ['0/donor.ccp4', '0/apolar.ccp4', '0/acceptor.ccp4']
grids = ['donor', 'apolar', 'acceptor']
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
cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(102.0), float(102.0), float(78.5), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([102.0,102.0,78.5], [105.060856925,103.12218675,76.5282538422], color="blue red", name="Arrows_rank_1_1")

cluster_dict["rank_1"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(100.655253837), float(100.131050767), float(79.3789846517), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(100.5), float(98.5), float(82.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([100.5,98.5,82.0], [103.45387441,98.4033154843,84.5691754425], color="red blue", name="Arrows_rank_1_2")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(104.0), float(101.0), float(79.5), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([104.0,101.0,79.5], [107.561829715,103.44217802,77.543243946], color="red blue", name="Arrows_rank_1_3")

cmd.load_cgo(cluster_dict["rank_1"], "Features_rank_1", 1)
cmd.load_cgo(cluster_dict["rank_1_arrows"], "Arrows_rank_1")
cmd.set("transparency", 0.2,"Features_rank_1")
cmd.group("Pharmacophore_rank_1", members="Features_rank_1")
cmd.group("Pharmacophore_rank_1", members="Arrows_rank_1")

if dirpath:
    f = join(dirpath, "0\label_threshold_rank_1.mol2")
else:
    f = "0\label_threshold_rank_1.mol2"

cmd.load(f, 'label_threshold_rank_1')
cmd.hide('everything', 'label_threshold_rank_1')
cmd.label("label_threshold_rank_1", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_rank_1', members= 'label_threshold_rank_1')


if dirpath:
    f = join(dirpath, '0\mesh.grd')
else:
    f = '0\mesh.grd'
cmd.load(f, 'mesh_0')
cmd.isomesh("isomesh_0", "mesh_0", 0.9)
cmd.color("grey80", "isomesh_0")
cmd.set('transparency', 0.4, "isomesh_0")

cmd.group('hotspot_0', "isomesh_0")
cmd.group('hotspot_0', "mesh_0")

if dirpath:
    f = join(dirpath, "1\label_threshold_15.2.mol2")
else:
    f = "1\label_threshold_15.2.mol2"

cmd.load(f, 'label_threshold_15.2')
cmd.hide('everything', 'label_threshold_15.2')
cmd.label("label_threshold_15.2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [15.2]
gfiles = ['1/donor.ccp4', '1/apolar.ccp4', '1/acceptor.ccp4']
grids = ['donor', 'apolar', 'acceptor']
num = 1
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

cluster_dict = {"rank_2":[], "rank_2_arrows":[]}
cluster_dict["rank_2"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(92.1223551058), float(73.7483900644), float(67.4622815087), float(1.0)]


cmd.load_cgo(cluster_dict["rank_2"], "Features_rank_2", 1)
cmd.load_cgo(cluster_dict["rank_2_arrows"], "Arrows_rank_2")
cmd.set("transparency", 0.2,"Features_rank_2")
cmd.group("Pharmacophore_rank_2", members="Features_rank_2")
cmd.group("Pharmacophore_rank_2", members="Arrows_rank_2")

if dirpath:
    f = join(dirpath, "1\label_threshold_rank_2.mol2")
else:
    f = "1\label_threshold_rank_2.mol2"

cmd.load(f, 'label_threshold_rank_2')
cmd.hide('everything', 'label_threshold_rank_2')
cmd.label("label_threshold_rank_2", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_rank_2', members= 'label_threshold_rank_2')


if dirpath:
    f = join(dirpath, '1\mesh.grd')
else:
    f = '1\mesh.grd'
cmd.load(f, 'mesh_1')
cmd.isomesh("isomesh_1", "mesh_1", 0.9)
cmd.color("grey80", "isomesh_1")
cmd.set('transparency', 0.4, "isomesh_1")

cmd.group('hotspot_1', "isomesh_1")
cmd.group('hotspot_1', "mesh_1")

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
