from __future__ import print_function
import os

if not os.path.exists("./hotspot_guided_docking"):
    os.mkdir("./hotspot_guided_docking")
os.system("""python ./docking.py "./hotspot_guided_docking" "1hcl" "ATP" """)

if not os.path.exists("./default_docking"):
    os.mkdir("./default_docking")
os.system("""python ./docking.py "./default_docking" "1hcl" "ATP" -hs False """)

print("Hotspot Guided")
os.system("""python ./comparision.py "1hcl" "1B38" "ATP" "hotspot_guided_docking/results.mol2" """)

print("Default")
os.system("""python ./comparision.py "1hcl" "1B38" "ATP" "default_docking/results.mol2" """)
