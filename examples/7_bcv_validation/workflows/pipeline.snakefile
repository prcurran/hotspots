
rule bcv_validation_devel:
    input:
        expand('data/{ypdb}/job.sh', ypdb=map(ypdb, all_apos)) +
        expand(

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

    output:
        done = 'data/{pdbpref}/{pdbid}/done.txt'

    run:
        with open(input.jobs, 'r') as f:
            cmd = f.read()

        print(cmd)

        with open(output.done, 'w') as w:
            w.write(cmd)
