
rule bcv_validation_devel:
    input:
        jobs = expand('data/{ypdb}/job.sh', ypdb=map(ypdb, all_apos))

rule make_job:
    output:
        expand('data/{ypdb}/job.sh', ypdb=map(xpdb, all_apos))

    run:
        df = pd.read_csv("results/inputs.csv")
        ligands = list(df.loc[df['apo'] == '{wildcard.pdbid}']['fragment_ID']) + list(df.loc[df['apo'] == '{wildcard.pdbid}']['lead_ID'])
        proteins = list(df.loc[df['apo'] == '{wildcard.pdbid}']['fragment']) + list(df.loc[df['apo'] == '{wildcard.pdbid}']['lead'])
        cmd = "python pipeline.py {} {} {}".format(pdb, ",".join(proteins), ",".join(ligands))
        with open(output, 'w') as f:
            f.write(cmd)

rule submit:
    input:
        jobs = expand('data/{ypdb}', ypdb=map(ypdb, all_apos))

    run:
        print(open('data/{wildcard.pdbpref}/{pdbid}/job.sh').read())





