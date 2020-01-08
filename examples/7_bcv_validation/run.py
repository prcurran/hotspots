"""
A Chris-proof run script
"""

import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor


#########################
workers = 2
#########################


def job(cmd):
    os.system(cmd)


def main():
    df = pd.read_csv("results/inputs.csv")
    hot_pdbs = set(df['apo'])
    vars = []
    for i, pdb in enumerate(hot_pdbs):
        ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
        proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])
        vars.append("python pipeline.py {} {} {}".format(pdb, ",".join(proteins), ",".join(ligands)))
        #job("python pipeline.py {} {} {}".format(pdb, ",".join(proteins), ",".join(ligands)))
    print vars
    with ProcessPoolExecutor(max_workers=workers) as executor:
        executor.map(job, vars)


if __name__ == "__main__":
    main()
