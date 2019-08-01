import pandas as pd
import seaborn as sns
import matplotlib.pyplot as ply

plt.rcParams['figure.figsize']=(10,20)

df = pd.read_csv("analysis.csv")
print(df)

a = df.loc[df["bcv_calc"] == 1]
b = df.loc[df["top_cavity"] == True]

ax = sns.barplot(x=b["ligand_id"], y=b["volume_overlap"], data=b)

for item in ax.get_xticklabels():
    item.set_rotation(90)

fig = ax.get_figure()
fig.savefig("figure3.png")