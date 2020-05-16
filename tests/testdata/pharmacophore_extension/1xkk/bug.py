from ccdc.protein import Protein
from ccdc import io


p0 = Protein.from_file(("1xkk.pdb"))
p1 = p0.copy()
p1.remove_all_waters()

with io.MoleculeWriter("1xkk_dry.pdb") as w:
    w.write(p1)

with io.MoleculeWriter("1xkk_dry.mol2") as w:
    w.write(p1)

p2 = Protein.from_file("1xkk_dry.pdb")
p3 = Protein.from_file("1xkk_dry.mol2")


for p, l in zip([p0, p2, p3], ["original", "mol2writer", "pdbwriter"]):
    print(l)
    print(f"# atoms: {len(p.atoms)}")
    print(f"# residues: {len(p.residues)}, # waters: {len(p.waters)}")
    print(p.residues[30:50])
    print(" ")
