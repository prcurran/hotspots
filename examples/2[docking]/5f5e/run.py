from __future__ import print_function
import os


prot = "5f5e"
ligand = "SAH"

ref = "5f5e"


if not os.path.exists("./hotspot_guided_docking"):
    os.mkdir("./hotspot_guided_docking")
os.system("""python ./docking.py "./hotspot_guided_docking" "{}" "{}" """.format(prot, ligand))

if not os.path.exists("./default_docking"):
    os.mkdir("./default_docking")
os.system("""python ./docking.py "./default_docking" "{}" "{}" -hs False """.format(prot, ligand))

# print("Hotspot Guided")
# os.system("""python ./comparision.py "{}" "{}" "{}" "hotspot_guided_docking/results.mol2" """.format(prot, ref, ligand))
#
# print("Default")
# os.system("""python ./comparision.py "{}" "{}" "{}" "default_docking/results.mol2" """.format(prot, ref, ligand))
