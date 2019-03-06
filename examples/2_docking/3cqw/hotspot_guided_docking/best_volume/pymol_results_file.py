
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
gfiles = ['apolar.grd', 'acceptor.grd']
grids = ['apolar', 'acceptor']
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


cluster_dict = {"18.4072484054":[], "18.4072484054_arrows":[]}

cluster_dict["18.4072484054"] += [COLOR, 1.00, 1.000, 0.000] + [ALPHA, 0.6] + [SPHERE, float(10.9629543701), float(-5.86685732217), float(19.2633291236), float(1.0)]


cluster_dict["18.4072484054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(9.5), float(-8.5), float(15.5), float(1.0)]

cluster_dict["18.4072484054_arrows"] += cgo_arrow([9.5,-8.5,15.5], [9.42807094925,-11.3091196945,14.7680681012], color="red blue", name="Arrows_18.4072484054_1")

cluster_dict["18.4072484054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(10.5), float(-3.0), float(24.0), float(1.0)]

cluster_dict["18.4072484054_arrows"] += cgo_arrow([10.5,-3.0,24.0], [11.4022665844,-4.07904317214,25.9751197812], color="red blue", name="Arrows_18.4072484054_2")

cluster_dict["18.4072484054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(11.5), float(-5.0), float(17.0), float(1.0)]

cluster_dict["18.4072484054_arrows"] += cgo_arrow([11.5,-5.0,17.0], [7.46517607067,-3.64803861043,17.0750787397], color="red blue", name="Arrows_18.4072484054_3")

cluster_dict["18.4072484054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-8.5), float(24.0), float(1.0)]

cluster_dict["18.4072484054_arrows"] += cgo_arrow([12.5,-8.5,24.0], [13.9912192931,-7.63608081942,27.03312466], color="red blue", name="Arrows_18.4072484054_4")

cluster_dict["18.4072484054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-0.5), float(15.0), float(1.0)]

cluster_dict["18.4072484054_arrows"] += cgo_arrow([12.5,-0.5,15.0], [13.2989121801,-1.5810167333,12.3970571675], color="red blue", name="Arrows_18.4072484054_5")

cluster_dict["18.4072484054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(12.5), float(-6.0), float(21.5), float(1.0)]

cluster_dict["18.4072484054_arrows"] += cgo_arrow([12.5,-6.0,21.5], [11.4022665844,-4.07904317214,25.9751197812], color="red blue", name="Arrows_18.4072484054_6")

cluster_dict["18.4072484054"] += [COLOR, 1.00, 0.00, 0.00] + [ALPHA, 0.6] + [SPHERE, float(14.0), float(-2.0), float(24.0), float(1.0)]

cluster_dict["18.4072484054_arrows"] += cgo_arrow([14.0,-2.0,24.0], [14.1992019713,-4.12704368017,26.5111222529], color="red blue", name="Arrows_18.4072484054_7")

cmd.load_cgo(cluster_dict["18.4072484054"], "Features_18.4072484054", 1)
cmd.load_cgo(cluster_dict["18.4072484054_arrows"], "Arrows_18.4072484054")
cmd.set("transparency", 0.2,"Features_18.4072484054")
cmd.group("Pharmacophore_18.4072484054", members="Features_18.4072484054")
cmd.group("Pharmacophore_18.4072484054", members="Arrows_18.4072484054")

if dirpath:
    f = join(dirpath, "label_threshold_18.4072484054.mol2")
else:
    f = "label_threshold_18.4072484054.mol2"

cmd.load(f, 'label_threshold_18.4072484054')
cmd.hide('everything', 'label_threshold_18.4072484054')
cmd.label("label_threshold_18.4072484054", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


cmd.group('Pharmacophore_18.4072484054', members= 'label_threshold_18.4072484054')

cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
