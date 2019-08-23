
rule bcv_validation_devel:
    input:
        expand('data/{ypdb}/job.sh', ypdb=map(ypdb, all_apos))

rule make_job:
    output:
        out_file = 'data/{pdbpref}/{pdbid}/job.sh'

    run:
        df = pd.read_csv("results/inputs.csv")
        ligands = list(df.loc[df['apo'] == wildcards.pdbid]['fragment_ID']) + list(df.loc[df['apo'] == wildcards.pdbid]['lead_ID'])
        proteins = list(df.loc[df['apo'] == wildcards.pdbid]['fragment']) + list(df.loc[df['apo'] == wildcards.pdbid]['lead'])
        cmd = "python pipeline.py {} {} {}".format(wildcards.pdbid, ",".join(proteins), ",".join(ligands))
        with open(output.out_file, 'w') as f:
            f.write(cmd)

rule submit:
    input:
        jobs = 'data/{pdbpref}/{pdbid}/job.sh'

    run:
        print(open('data/{wildcards.pdbpref}/{wildcards.pdbid}/job.sh').read())





