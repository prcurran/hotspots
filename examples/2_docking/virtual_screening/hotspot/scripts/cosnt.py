from ccdc import io
from hotspots.hs_io import HotspotReader


hs = HotspotReader("out.zip").read()
c = hs.docking_constraint_atoms(max_constraints=1)
d = c.to_molecule()

with io.MoleculeWriter("constraints.mol2") as w:
    w.write(d)