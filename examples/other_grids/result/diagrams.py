from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw



from ccdc.io import EntryReader
import pandas as pd


mols = "./pharmit/query_results.sdf"
mols = EntryReader(mols)


smiles = []
name = []
rmsd = []

for m in mols:
    smiles.append(m.molecule.smiles)
    name.append(m.identifier)
    rmsd.append(m.attributes["rmsd"])

df = pd.DataFrame({"smiles":smiles, "name":name, "rmsd":rmsd})

df.to_csv("pharmit.csv")

# entries= [(m.smiles, m.identifier) for m in mols]
#
# ligs =[]
#
# for i, entry in enumerate(entries):
#     lig = Chem.MolFromSmiles(entry[0])
#     name = ""
#     for e in entry[1].split(" "):
#         name += "{}\n".format(e)
#
#     lig.SetProp("_Name", "{} \n {}".format(str(i), name))
#     ligs.append(lig)
#
# for i, l in enumerate(ligs):
#     tmp = AllChem.Compute2DCoords(l)
#
# img = Draw.MolsToGridImage(mols=ligs[:8],
#                      molsPerRow=4,
#                      subImgSize=(200, 200),
#                      legends=[l.GetProp("_Name") for l in ligs[:8]])
#
# img.save("mols.png")