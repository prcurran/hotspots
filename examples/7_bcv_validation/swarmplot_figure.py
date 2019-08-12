import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


df = pd.read_csv("analysis.csv")

df = df.loc[df["ligand_cavity"] == True]
df = df.loc[df["atom_type"] != "apolar"]

buriedness_method = ["ligsite", "ghecom", "ghecom_internal"]

pal = sns.color_palette(["#3380FF", "#FFC433", "#FC1B22"])

ax = sns.swarmplot(x=df["buriedness_method"], y=df["atomic_overlap"], hue=df["atom_type"],
                   hue_order=["donor", "apolar", "acceptor"], data=df, palette=pal)

ax.legend_.remove()

ax.set(xlabel='Buriedness Method', ylabel="Percentage Atomic Volume Overlap", title="", ylim=[-5, 105])

fig = ax.get_figure()
fig.savefig("swarmplot1_polar.png")