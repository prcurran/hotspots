
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

cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(101.5), float(100.5), float(79.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([101.5,100.5,79.0], [99.3539190176,97.9993265054,75.8282606673], color="blue red", name="Arrows_rank_1_1")

cluster_dict["rank_1"] += [COLOR, 0.00, 0.00, 1.00] + [ALPHA, 0.6] + [SPHERE, float(101.5), float(100.5), float(79.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([101.5,100.5,79.0], [99.3539190176,97.9993265054,75.8282606673], color="blue red", name="Arrows_rank_1_2")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(100.0), float(93.0), float(83.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([100.0,93.0,83.0], [102.239887618,91.3185087631,82.099199525], color="red blue", name="Arrows_rank_1_3")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(100.0), float(103.0), float(76.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([100.0,103.0,76.0], [99.3359192134,100.912247039,74.459274015], color="red blue", name="Arrows_rank_1_4")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(100.5), float(100.5), float(77.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([100.5,100.5,77.0], [99.3359192134,100.912247039,74.459274015], color="red blue", name="Arrows_rank_1_5")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(101.0), float(98.5), float(83.0), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([101.0,98.5,83.0], [103.45387441,98.4033154843,84.5691754425], color="red blue", name="Arrows_rank_1_6")

cluster_dict["rank_1"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(105.5), float(102.5), float(79.5), float(1.0)]

cluster_dict["rank_1_arrows"] += cgo_arrow([105.5,102.5,79.5], [107.561829715,103.44217802,77.543243946], color="red blue", name="Arrows_rank_1_7")

cluster_dict["rank_1"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(96.4404638051), float(92.8539517401), float(82.3560117005), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(96.5), float(93.0), float(84.0), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(101.318937137), float(99.7879379828), float(80.3960782506), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(97.0), float(96.5), float(80.0), float(1.0)]


cluster_dict["rank_1"] += [COLOR, 0.6, 0.1, 0.6] + [ALPHA, 0.6] + [SPHERE, float(97.5), float(99.5), float(76.0), float(1.0)]


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
