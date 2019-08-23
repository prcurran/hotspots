
rule bcv_validation_devel:
    input:
        df = "results/inputs.csv"
        apo = expand('data/{xpdb}/{pdbid}.pdb', xpdb=map(xpdb, all_pdb_codes))

rule make_job:
    output:
        expand('data/{xpdb}/job.sh', xpdb=map(xpdb, all_pdb_codes))

    run:
        df = pd.read_csv(input.df)
        ligands = list(df.loc[df['apo'] == '{wildcard.pdbid}']['fragment_ID']) +
                  list(df.loc[df['apo'] == '{wildcard.pdbid}']['lead_ID'])

        proteins = list(df.loc[df['apo'] == '{wildcard.pdbid}']['fragment']) +
                   list(df.loc[df['apo'] == '{wildcard.pdbid}']['lead'])
        cmd = "python pipeline.py {} {} {}".format(pdb, ",".join(proteins), ",".join(ligands))
        with open(output, 'w') as f:
            f.write(cmd)
