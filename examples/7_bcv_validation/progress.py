"""

track the progress of the BCV validation

"""
from __future__ import print_function
import datetime

import pandas as pd

from pipeline import HotspotPipeline


prefix = "/vagrant/github_pkgs/hotspots/examples/7_bcv_validation"
buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']


def main():
    df = pd.read_csv("inputs.csv")
    hot_pdbs = set(df['apo'])
    reports = []

    for i, pdb in enumerate(hot_pdbs):
        for method in buriedness_methods:
            ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
            proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])

            hp = HotspotPipeline(apo=pdb, buriedness_method=method, protein_id=proteins, ligand_id=ligands)
            report = hp.status()
            reports.append(report)

    progress = pd.concat(reports)

    date = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    progress.to_csv("progress_report_{}.csv".format(date))
    print(progress)


if __name__ == "__main__":
    main()