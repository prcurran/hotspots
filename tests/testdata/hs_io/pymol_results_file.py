
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

dirpath = None
    
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

cmd.load(f, 'label_threshold_10_0')
cmd.hide('everything', 'label_threshold_10_0')
cmd.label("label_threshold_10_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


if dirpath:
    f = join(dirpath, "0/label_threshold_14_0.mol2")
else:
    f = "0/label_threshold_14_0.mol2"

cmd.load(f, 'label_threshold_14_0')
cmd.hide('everything', 'label_threshold_14_0')
cmd.label("label_threshold_14_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


if dirpath:
    f = join(dirpath, "0/label_threshold_17_0.mol2")
else:
    f = "0/label_threshold_17_0.mol2"

cmd.load(f, 'label_threshold_17_0')
cmd.hide('everything', 'label_threshold_17_0')
cmd.label("label_threshold_17_0", "name")
cmd.set("label_font_id", 7)
cmd.set("label_size", -0.4)


colour_dict = {'acceptor':'red', 'donor':'blue', 'apolar':'yellow', 'negative':'purple', 'positive':'cyan'}

threshold_list = [10, 14, 17]
gfiles = ['0/acceptor_0.grd', '0/apolar_0.grd', '0/donor_0.grd']
grids = ['acceptor', 'apolar', 'donor']
num = '0'
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
            cmd.group('threshold_%s_%s' % (t, num), members = 'surface_%s_%s_%s'%(grids[i],t, num))
            cmd.group('threshold_%s_%s' % (t, num), members='label_threshold_%s_%s' % (t, num))
        except:
            continue



    try:
        cmd.group('hotspot_%s' % (num), members='threshold_%s_%s' % (t, num))
    except:
        continue
    
    for g in grids:
        
        cmd.group('hotspot_%s' % (num), members='%s_%s' % (g,num))


cmd.bg_color("white")
cmd.show("cartoon", "protein")
cmd.color("slate", "protein")
cmd.show("sticks", "organic")
cmd.hide("lines", "protein")
