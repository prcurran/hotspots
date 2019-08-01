"""

generate paper figures

"""

# buriedness_method, cavity_id, prot_id, lig_id, top (bool), volume_overlap


from __future__ import print_function
import datetime

import pandas as pd

from pipeline import HotspotPipeline


prefix = "/vagrant/github_pkgs/hotspots/examples/7_bcv_validation"
buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']


def main():
    df = pd.read_csv("inputs.csv")
    frags = set(df['fragment'])
    leads = set(df['lead'])

    hot_pdbs = set(df['apo'])
    reports = []

    for i, pdb in enumerate(hot_pdbs):
        for method in buriedness_methods:
            ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
            proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])

            hp = HotspotPipeline(apo=pdb, buriedness_method=method, protein_id=proteins, ligand_id=ligands)
            report = hp.analysis()
            reports.append(report)
            print(report)

    dat = pd.concat(reports, ignore_index=True)

    classification = []
    for a in list(dat['other_id']):
        if a in frags:
            classification.append("fragment")
        elif a in leads:
            classification.append("lead")

    dat["lig_class"] = classification

    dat.to_csv("analysis.csv")


if __name__ == "__main__":
    main()