from hotspots.hs_pharmacophore import PharmacophoreModel
from hotspots import hs_io

from hotspots.result import Results

import os
from ccdc import io


mode = "ligand_pose_comparison"
base = "/home/pcurran/question1"

targets = {
           #  "CDK2": ["1aq1"],
           # "DHFR": ["1drf"],
           # "Thrombin": ["1c4v"],
           # "HIVRT": ["1tvr"],
           "A2Ar": ["2ydo"]}

chains = {"1hcl": "A",
          "1aq1": "A",
          "1drf": "A",
          "2w9t": "A",
          "1c4v": "2",
          "1vr1": "H",
          "1tvr": "A",
          "1dlo": "A",
          "2ydo": "A",
          }

ligands = {"1aq1": "STU_A_299",
          "1drf": "FOL_A_187",
          "1c4v": "IH2_2_370",
          "1tvr": "TB9_A_600",
          "2ydo": "ADN_A_400",
          }

for target, pdbs in targets.items():
    print target
    for pdb in pdbs:
        chain = chains[pdb]
        ligand_id = ligands[pdb]

        out_dir = os.path.join(base, target, pdb, "reference")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        try:
            p = PharmacophoreModel._from_siena(pdb, ligand_id, mode, target, out_dir=out_dir)
            p.write(os.path.join(out_dir, "reference_pharmacophore.py"))


            prot = hs_io.HotspotReader(os.path.join(base, target, pdb,"out.zip")).read().protein

            hs = Results(protein=prot, super_grids=p.dic)

            with hs_io.HotspotWriter(out_dir) as wf:
                wf.write(hs)

            with io.MoleculeWriter(os.path.join(out_dir, "aligned.mol2")) as w:
                for l in p.representatives:
                    w.write(l)
        except RuntimeError:
            print "skipped {}".format(target)
