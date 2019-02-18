
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
gfiles = ['donor.grd', 'apolar.grd', 'acceptor.grd']
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


cluster_dict = {"17.7479991913":[], "17.7479991913_arrows":[]}

cluster_dict["17.7479991913"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.5), float(-1.5), float(53.0), float(1.0)]

cluster_dict["17.7479991913_arrows"] += cgo_arrow([28.5,-1.5,53.0], [27.4034055644,-3.1090460132,55.890827172], color="blue red", name="Arrows_17.7479991913_1")

cluster_dict["17.7479991913"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(28.0), float(2.5), float(54.0), float(1.0)]

cluster_dict["17.7479991913_arrows"] += cgo_arrow([28.0,2.5,54.0], [27.5574078436,3.0720454656,56.8888419424], color="blue red", name="Arrows_17.7479991913_2")

cluster_dict["17.7479991913"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(24.766061346), float(0.580047837555), float(52.8181029438), float(1.0)]


cluster_dict["17.7479991913"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(22.0), float(-2.0), float(55.0), float(1.0)]

cluster_dict["17.7479991913_arrows"] += cgo_arrow([22.0,-2.0,55.0], [20.8983092904,-2.8660424168,57.745854626], color="red blue", name="Arrows_17.7479991913_3")

cluster_dict["17.7479991913"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(23.0), float(1.5), float(50.5), float(1.0)]


cluster_dict["17.7479991913"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(24.5), float(4.0), float(53.5), float(1.0)]

cluster_dict["17.7479991913_arrows"] += cgo_arrow([24.5,4.0,53.5], [26.8273970396,5.510081548,53.2487880704], color="red blue", name="Arrows_17.7479991913_4")

cluster_dict["17.7479991913"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(1.5), float(53.5), float(1.0)]

cluster_dict["17.7479991913_arrows"] += cgo_arrow([29.5,1.5,53.5], [32.0534743844,2.620038776,52.6737795604], color="red blue", name="Arrows_17.7479991913_5")

cluster_dict["17.7479991913"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(29.5), float(1.5), float(53.5), float(1.0)]

cluster_dict["17.7479991913_arrows"] += cgo_arrow([29.5,1.5,53.5], [32.0534743844,2.620038776,52.6737795604], color="red blue", name="Arrows_17.7479991913_6")

cmd.load_cgo(cluster_dict["17.7479991913"], "Features_17.7479991913", 1)
cmd.load_cgo(cluster_dict["17.7479991913_arrows"], "Arrows_17.7479991913")
cmd.set("transparency", 0.2,"Features_17.7479991913")
cmd.group("Pharmacophore_17.7479991913", members="Features_17.7479991913")
cmd.group("Pharmacophore_17.7479991913", members="Arrows_17.7479991913")

if dirpath:
    f = join(dirpath, "label_threshold_17.7479991913.mol2")
else:
    f = "label_threshold_17.7479991913.mol2"

cmd.load(f, 'label_threshold_17.7479991913')
cmd.hide('everything', 'label_threshold_17.7479991913')
cmd.label("label_threshold_17.7479991913", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_17.7479991913', members= 'label_threshold_17.7479991913')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
