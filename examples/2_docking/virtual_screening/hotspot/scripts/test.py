from hotspots import calculation
from hotspots.hs_io import HotspotWriter

r = calculation.Runner()

result = r.from_pdb("3cqw", charged_probes=False, nprocesses=3)
with HotspotWriter("/home/pcurran/New folder/akt1/protoss") as w:
    w.write(result)

print set([len(a.neighbours) for a in result.protein.atoms])