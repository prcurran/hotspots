import os
import pandas as pd


def all_pdbs():
    f = "results/inputs.csv"
    df = pd.read_csv(f)

    return list(set(df['apo']))

def xpdb(pdbid): return '{0}/{1}/{1}.pdb'.format(pdbpref(pdbid), pdbid) # e.g. '2fdu' => 'fd/2fdu'

def ypdb(pdbid): return '{0}/{1}'.format(pdbpref(pdbid), pdbid) # e.g. '2fdu' => 'fd/2fdu'

def pdbpref(pdbid): return pdbid[1:3] # e.g. '2fdu' => 'fd'

def pdbbase(pdbid): return '%s/pdb%s' % (pdbpref(pdbid), pdbid) # e.g. '2fdu' => 'fd/2fdu'

def fpdb(pdb_id, step, suffix):
    if pdb_id is None:
        pdb_id = '{pdbid}'
        pdb_pref = '{pdbpref}'
    else:
        pdb_pref = pdbpref(pdb_id)
    return '%s/%s/%s%s' % (step, pdb_pref, pdb_id, suffix)

all_apos = all_pdbs()