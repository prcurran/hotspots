import pandas as pd
import numpy as np
import operator
import seaborn as sns; sns.set(style="white", color_codes=True)
import matplotlib.pyplot as plt
import os
from ccdc.descriptors import StatisticalDescriptors


def process_data(f):
    lines = [l.strip("\n") for l in open(f, "rb").readlines()]

    header = lines[5]
    header = [b.strip() for b in [a for a in header.split("  ") if a != '' and a != '#']]
    print header

    all = lines[7:]
    cat = zip(*[[a for a in entry.split(" ") if a != ''] for entry in all])
    df = pd.DataFrame({h: cat[i] for i, h in enumerate(header)})

    df["actives"] = np.array(map(lambda x: 'CHEMBL' in x, list(df['Ligand name'])))
    return df.sort_values(by=['Score'], ascending=False)

#########################################################################
base = "/home/pcurran/github_packages/hotspots/examples/2_docking/virtual_screening/search_efficiency_100/akt1"
#########################################################################

runs = [0, 10, 100]
erich_levels = [0.01,0.05,0.1]
frames = []

weight = []
ef = []
efl = []

for r in runs:
    f = os.path.join(base, "docking/{}/bestranking.lst".format(r))
    df = process_data(f)

    new = zip(list(df["Ligand name"]), list(df["Score"]), list(df["actives"]))

    sd = StatisticalDescriptors()

    rank = sd.RankStatistics(new, operator.itemgetter(2))
    tpr, fpr = rank.ROC()

    efl.extend(erich_levels)
    ef.extend([rank.EF(e) for e in erich_levels])
    weight.extend([r]*len(erich_levels))

    df['tpr'] = tpr
    df['fpr'] = fpr
    df['constraint'] = [int(r)] * len(tpr)
    if int(r) == 0:
        df['DE(con'] = [0.0] * len(tpr)

    frames.append(df)

con = pd.concat(frames, sort=True)
print con
palette = sns.color_palette("Reds", 3)
ax = sns.lineplot(x=con['fpr'], y=con['tpr'], hue=con['constraint'], palette=palette)
sns.lineplot(x=[0,1], y=[0,1], color="grey")
plt.title("Docking Enrichment for AKT1 (PDB:'3cqw')")
plt.savefig(os.path.join(base, "roc.png"))
plt.close()

# timer = pd.DataFrame({'constraint': np.array(con['constraint']).astype(int),
#                       'time': np.array(con['time']).astype(int)}, index=list(con['Ligand name']))

con['time'] = [float(x) for x in con['time']]
sns.barplot(x=con['constraint'], y=con['time'], palette=palette)
plt.title("Runtime(seconds) per Ligand for AKT1 (PDB:'3cqw')")
plt.savefig(os.path.join(base, "Time.png"))
plt.close()

enr_df = pd.DataFrame({"Weight": weight, "Enrichment Factor": ef, "EF Level": efl})
ax = sns.barplot(x=enr_df['EF Level'], y=enr_df['Enrichment Factor'], hue= enr_df['Weight'], palette=palette)
plt.title("Enrichment Factors for AKT1 (PDB:'3cqw')")
plt.savefig(os.path.join(base,"EF.png"))
plt.close()


