import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def _add_classification(df):

    i = pd.read_csv("inputs.csv")
    frags = set(i['fragment'])
    leads = set(i['lead'])
    classification = []

    for a in list(df['other_id']):
        if a in frags:
            classification.append("fragment")
        elif a in leads:
            classification.append("lead")

    df["classification"] = classification
    return df


def barplot1(df):

    ligs = set(df["apo"])
    order_dic = {np.mean([float(a) for a in list(df.loc[df["apo"] == l]["volume_overlap"])]): l for l in ligs}

    ks = sorted(order_dic.keys(), reverse=True)
    ord = [order_dic[k] for k in ks]

    g = sns.FacetGrid(df, col="buriedness_method", height=7, aspect=1, legend_out=True)
    g = (g.map(sns.barplot, 'apo', 'volume_overlap', 'lig_class', order=ord, hue_order=["fragment", "lead"],
               ci=None, palette="Reds").add_legend())

    g.set_xticklabels(rotation=45)
    g.set(xlabel='', ylabel="", title="")

    plt.savefig("barplot1.png")
    plt.close()


def barplot2(df, clas = "fragment"):
    def over_threshold(l, threshold):
        return [l for item in l if item >= threshold]

    df = df.loc[df["atom_type"] != 'apolar']
    df = df.loc[df["classification"] == clas]

    buriedness_method = ["ligsite", "ghecom", "ghecom_internal"]
    data = {}

    for b in buriedness_method:
        x = df.loc[df["buriedness_method"] == b]
        x = x.loc[x["atom_type"] != "apolar"]
        data.update({b: list(x["atomic_overlap"])})

    total = max([len(k) for k in data.keys()])

    bm = []
    perc = []
    thres = []

    thresholds = [1, 5, 10, 50, 100]
    for t in thresholds:
        for b, d in data.items():
            bm.append(b)
            thres.append(t)
            perc.append((len(over_threshold(d, t)) / total) * 100)

    ndf = pd.DataFrame({"thresholds": thres, "passed": perc, "buriedness_method": bm})

    ax = sns.barplot(x=ndf["thresholds"], y=ndf["passed"], hue=ndf["buriedness_method"],
                     hue_order= ['ligsite', 'ghecom', 'ghecom_internal'], palette='Reds')

    ax.set(xlabel="Overlap threshold (%)", ylabel="Atomic overlap greater than threshold (%)", title="", ylim=[0, 100])
    plt.savefig("barplot2_{}.png".format(clas))
    plt.close()


def boxplot1(df):

    ax = sns.boxplot(x='buriedness_method', y='volume_overlap', hue="classification", data=df, order=["ligsite", "ghecom", "ghecom_internal"],
                palette="Reds")

    ax.set(xlabel='Buriedness Method', ylabel="Percentage Volume Overlap", title="", ylim=[-5,100])

    plt.savefig("boxplot1.png")
    plt.close()


if __name__ == "__main__":
    sns.set(style="ticks", color_codes=True)
    df = pd.read_csv("analysis.csv")
    df = df.loc[df["ligand_cavity"] == True]
    # barplot1(df)
    # barplot1(df)
    barplot2(df, clas="fragment")
    barplot2(df, clas="lead")
