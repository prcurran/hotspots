from __future__ import division
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Initialise seaborn styling
sns.set(style="white",
        rc={"axes.facecolor": (0, 0, 0, 0)},
        font_scale=6)
# pal = [(117/255, 216/255, 213/255),
#        (211/255, 60/255, 96/255)]

pal = [(202/255, 41/255, 54/255),
       (140/255, 157/255, 196/255)]

# Get Data
fname = "scores.csv"
df = pd.read_csv(fname)
df["median"] = 0
pdbs = set(df['pdb'])
for p in pdbs:
    median = np.median(df.loc[df['pdb'] == p]["scores"])
    a = df.loc[df['pdb'] == p, 'median'] = median
df = df.sort_values(by='median', ascending=False)

# Initialize the FacetGrid object
g = sns.FacetGrid(df, row="pdb", hue="tractability", height=4, aspect=75, size=.5, palette=pal)

# Draw the densities in a few steps
g.map(sns.kdeplot, "scores", clip_on=False, shade=True, alpha=.7, lw=3, bw=.2)
g.map(sns.kdeplot, "scores", clip_on=False, color="w", lw=9, bw=.2)
g.set(xlim=(0, 25))


# Set the subplots to overlap
g.fig.subplots_adjust(hspace=-.9)

# Remove axes details that don't play well with overlap
g.set_titles("")
g.set(yticks=[])
g.despine(bottom=True, left=True)

# Create legend, and position
t = {"d": "Druggable", "n": "Less-Druggable"}
labels = [t[s] for s in set(df["tractability"])]
colors = pal
handles= [patches.Patch(color=col, label=lab) for col, lab in zip(colors, labels)]
l = plt.legend(handles=handles, title='Tractability', loc="upper right", bbox_to_anchor=(0.3, 7.5))
frame = l.get_frame()
frame.set_facecolor('white')
frame.set_edgecolor('white')

plt.savefig("tractability.png")
