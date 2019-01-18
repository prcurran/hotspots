#!/bin/sh
mkdir ./hotspot_guided_docking
python ./docking.py "hotspot_guided_docking" "1hcl" "ATP" -hs "True"

mkdir ./default_docking
python ./docking.py "default_docking" "1hcl" "ATP" -hs "False"

echo "Hotspot Guided"
python ./comparision.py "1hcl" "1B38" "ATP" "hotspot_guided_docking/results.mol2"

echo "Default"
python ./comparision.py "1hcl" "1B38" "ATP" "default_docking/results.mol2"