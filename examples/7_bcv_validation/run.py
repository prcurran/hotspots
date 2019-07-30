"""
A Chris-proof run script
"""

import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


#########################
workers = 1
#########################


def job(cmd):
    os.system(cmd)


def main():
    df = pd.read_csv("inputs.csv")
    hot_pdbs = set(df['apo'])
    vars = []
    for i, pdb in enumerate(hot_pdbs):
        ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
        proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])
        vars.append("python pipeline.py {} {} {}".format(pdb, ",".join(proteins), ",".join(ligands)))

    with ThreadPoolExecutor(max_workers=workers) as executor:
        executor.submit(job, vars[0])


if __name__ == "__main__":
    main()
