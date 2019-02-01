import os

from ccdc.protein import Protein

from pdb_python_api import PDBResult
from hotspots.hs_pharmacophore_extension import PharmacophoreModel
from hotspots.calculation import Results
from hotspots.hs_io import HotspotWriter, HotspotReader


dirname = "./result"
pdb = "1vr1"
reps = "representatives.dat"

if not os.path.exists(dirname):
    os.mkdir(dirname)

PDBResult(identifier=pdb).download(out_dir=dirname)

if os.path.exists(reps):
    representatives = reps
else:
    representatives = None

try:
    result = HotspotReader(path=os.path.join(dirname, "out.zip")).read()
    pharmacophore = result.get_pharmacophore_model()
    pharmacophore.rank_features(max_features=5)

except:
    pharmacophore = PharmacophoreModel.from_pdb(pdb_code=pdb, chain="H", out_dir=dirname, representatives=representatives)
    pharmacophore.rank_features(max_features=5)
    result = Results(super_grids=pharmacophore.dic, protein=Protein.from_file(os.path.join(dirname, pdb + ".pdb")))


pharmacophore.write(os.path.join(dirname, "crossminer.cm"))
pharmacophore.write(os.path.join(dirname, "pharmit.json"))
# write out Results object
settings = HotspotWriter.Settings()
settings.isosurface_threshold=[2,5,10]
with HotspotWriter(dirname, settings=settings) as w:
    w.write(result)