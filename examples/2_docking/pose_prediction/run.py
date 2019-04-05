from __future__ import print_function
import os


#####################################################################
# input data
prot = "1m17"       # protein to dock into
ligand = "AQ4"      # ligand to dock
ref = "1m17"        # reference structure for comparison
#####################################################################

# run docking
# hotspot-guided GOLD
if not os.path.exists("./hotspot_guided_docking"):
    os.mkdir("./hotspot_guided_docking")

os.system("""python ./docking.py "./hotspot_guided_docking" "{}" "{}" """.format(prot, ligand))

# default GOLD
if not os.path.exists("./default_docking"):
    os.mkdir("./default_docking")

os.system("""python ./docking.py "./default_docking" "{}" "{}" -hs False """.format(prot, ligand))

# comparision
os.system("""python ./comparision.py "{}" "{}" "{}" "hotspot_guided_docking/results.mol2" """.format(prot, ref, ligand))
os.system("""python ./comparision.py "{}" "{}" "{}" "default_docking/results.mol2" """.format(prot, ref, ligand))
