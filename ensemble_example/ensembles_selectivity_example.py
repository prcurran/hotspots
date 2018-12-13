from __future__ import print_function, division
import os
from os.path import join, exists
from glob import glob
from hotspots.calculation import Results
from hotspots import hs_io



def make_savedir(prot_name, tag=None):
    if tag:
        save_dir = join(os.getcwd(), "{}_{}".format(tag, prot_name))
    else:
        save_dir = join(os.getcwd(),  prot_name)
    if not exists(save_dir):
        os.mkdir(save_dir)
    return save_dir
 
 
prot1_name= "BRD1"
prot2_name ="BAZ2B"

prot1_paths = glob(join(os.getcwd(),"{}*".format( prot1_name), "out.zip"))
print(prot1_paths)
prot2_paths = glob(join(os.getcwd(),"{}*".format( prot2_name), "out.zip"))
print(prot2_paths)

prot1_res_list = [hs_io.HotspotReader(p).read() for p in prot1_paths]
prot2_res_list = [hs_io.HotspotReader(p).read() for p in prot2_paths]

# Calculate ensemble hotspots for the two proteins
ensemble_1 = Results.from_grid_ensembles(prot1_res_list, prot1_name)
#Save ensemble:
out1 = make_savedir(prot1_name, "ensemble")
with hs_io.HotspotWriter(out1, visualisation="pymol", grid_extension=".ccp4", zip_results=True) as writer:
     writer.write(ensemble_1)
del(prot1_res_list)

ensemble_2 = Results.from_grid_ensembles(prot2_res_list, prot2_name)
#Save ensemble:
out2 = make_savedir(prot2_name, "ensemble")
with hs_io.HotspotWriter(out2, visualisation="pymol", grid_extension=".ccp4", zip_results=True) as writer:
     writer.write(ensemble_2)
del(prot2_res_list)

# Get selectivity maps
sel_map_result = ensemble_1.get_selectivity_map(ensemble_2, 1)
out_sel = make_savedir("select_{}_{}".format(prot1_name, prot2_name))
with hs_io.HotspotWriter(out_sel, visualisation="pymol", grid_extension=".ccp4", zip_results=True) as writer:
     writer.write(sel_map_result)
