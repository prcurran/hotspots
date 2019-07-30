import os
import pandas as pd


def main():
    # prefix = "/hps/nobackup/research/chembl/hotspots/bcv_validation"
    prefix = "/vagrant/github_pkgs/hotspots/examples/7_pipeline"
    df = pd.read_csv("inputs.csv")
    hot_pdbs = set(df['apo'])

    for pdb in hot_pdbs:
        ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
        proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])

        cmd = "python test_pipeline.py {} {} {} {}".format(
                                                        prefix,
                                                        pdb,
                                                        ",".join(proteins),
                                                        ",".join(ligands),
                                                        )

        print(cmd)
        os.system(cmd)


if __name__ == "__main__":
    main()
