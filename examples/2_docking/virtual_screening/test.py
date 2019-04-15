from hotspots import hs_docking
from hotspots import hs_io

from ccdc.io import MoleculeWriter

hs = hs_io.HotspotReader("/home/pcurran/github_packages/hotspots/examples/2_docking/virtual_screening/akt1/hotspot/out.zip").read()

c =hs.docking_constraint_atoms(max_constraints=3)
mol = c.to_molecule()

with MoleculeWriter("c.mol2") as w:
    w.write(mol)