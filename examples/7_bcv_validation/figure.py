import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def over_threshold(l, threshold):
    return [l for item in l if item >= threshold]


plt.rcParams['figure.figsize']=(8,10)
df = pd.read_csv("analysis.csv")

df = df.loc[df["ligand_cavity"] == True]
df = df.loc[df["lig_class"] == "lead"]
# df = df.loc[df["atom_type"] != "apolar"]

buriedness_method = ["ligsite", "ghecom", "ghecom_internal"]
data = {}

for b in buriedness_method:
    x = df.loc[df["buriedness_method"] == b]
    x = x.loc[x["atom_type"] != "apolar"]
    data.update({b: list(x["atomic_overlap"])})

total = 83

bm =[]
perc = []
thres = []

thresholds = [1, 5, 10, 50, 100]
for t in thresholds:
    for b, d in data.items():
        bm.append(b)
        thres.append(t)
        perc.append((len(over_threshold(d, t)) / total) * 100)


ndf = {"thresholds": thres, "passed": perc, "buriedness_method":bm}
ax = sns.barplot(x=ndf["thresholds"], y=ndf["passed"], hue=ndf["buriedness_method"])
ax.set_ylabel("Percentage of fragment atoms with overlap greater than threshold")
ax.set_xlabel("Overlap threshold")

fig = ax.get_figure()
fig.savefig("figure8.png")

# print(df)
# # donor acceptor apolar
#
# pal = sns.color_palette(["#3380FF", "#FFC433", "#FC1B22"])
#
# print(pal)
# # s = df.sort_values(by="max_overlap", ascending=0)
# #
# #
# # # b = df.loc[df["bcv_calc"] == 1]
# # # # b = df.loc[df["top_cavity"] == True]
# #
# ax = sns.swarmplot(x=df["buriedness_method"], y=df["atomic_overlap"], hue=df["atom_type"], data=df, palette=pal)
# #
# for item in ax.get_xticklabels():
#     item.set_rotation(90)
#
# fig = ax.get_figure()
# fig.savefig("figure6.png")