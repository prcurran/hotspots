import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as patches
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import seaborn as sns
import concurrent

from hotspots.calculation import Runner
from hotspots.result import Extractor
from hotspots.grid_extension import Grid
from ccdc.descriptors import StatisticalDescriptors as sd
import operator


def tractability_workflow(protein, tag):
    """
    A very simple tractability workflow.

    :param str protein: PDB identification code
    :param str tag: Tractability tag: either 'druggable' or 'less-druggable'
    :return: `pandas.DataFrame`
    """
    # 1) calculate Fragment Hotspot Result
    runner = Runner()
    result = runner.from_pdb(pdb_code=protein,
                             nprocesses=1,
                             buriedness_method='ghecom')

    # 2) calculate Best Continuous Volume
    extractor = Extractor(hr=result)
    bcv_result = extractor.extract_volume(volume=500)

    # 3) find the median score
    grid = Grid.single_grid(bcv_result.super_grids)
    values = grid.grid_values(threshold=5)
    median = np.median(values)

    # 4) return the data
    return pd.DataFrame({'scores': values,
                         'pdb': [protein] * len(values),
                         'median': [median] * len(values),
                         'tractability': [tag] * len(values),
                         })


def joyplot(df, fname='joy.png'):
    """
    Visualises the Fragment Hotspot Maps score distributions as a series of stacked kernal density estimates ordered by
    median value.

    Adapted from the seaborn gallery:
        https://seaborn.pydata.org/examples/kde_ridgeplot.html

    :param `pandas.DataFrame` df: Fragment Hotspot scores data
    :return None
    """
    sns.set(style="white",
            rc={"axes.facecolor": (0, 0, 0, 0)},
            font_scale=6)

    palette = ["#5bd9a4",
               "#c75048",
               "#808080"
               ]

    # Initialize the FacetGrid object
    ax = sns.FacetGrid(df, row="pdb", hue="tractability", height=2, aspect=40, palette=palette)

    # Draw the densities in a few steps
    ax.map(sns.kdeplot, "scores", clip_on=False, shade=True, alpha=.7, lw=3, bw=.2)
    ax.map(sns.kdeplot, "scores", clip_on=False, color="w", lw=9, bw=.2)

    # Format the plots
    ax.set(xlim=(0, 25))  #
    ax.fig.subplots_adjust(hspace=-.9)
    ax.set_titles("")
    ax.set(yticks=[])
    ax.despine(bottom=True, left=True)

    # Create legend, and position
    tag = {"d": "Druggable", "n": "Less-Druggable", "unknown": "Unknown"}
    labels = [tag[s] for s in set(df["tractability"])]
    handles = [patches.Patch(color=col, label=lab) for col, lab in zip(palette, labels)]
    legend = plt.legend(handles=handles, title='Tractability', loc="upper right", bbox_to_anchor=(0.3, 7.5))
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')

    plt.savefig(fname)
    plt.close()


def rocplot(data, a_col=1, fname='roc.png'):
    """
    Create a ROC Curve using seaborn

    :param lists data: supply ranked data as list of list
    :return: None
    """
    rs = sd.RankStatistics(scores=data, activity_column=operator.itemgetter(a_col))
    tpr, fpr = rs.ROC()
    ax = sns.lineplot(x=fpr, y=tpr, estimator=None, color="#c75048")
    ax = sns.lineplot(x=[0,1], y=[0,1], color="grey")
    ax.set(xlabel='FPR', ylabel='TPR', title=f"ROC Curve (AUC: {rs.AUC():.2f})")

    plt.savefig(fname)
    plt.close()


def main():
    subset = {'1e9x': 'd', '1udt': 'd', '2bxr': 'd', '1r9o': 'd', '3d4s': 'd', '1k8q': 'd',
              '1xm6': 'd', '1rwq': 'd', '1yvf': 'd', '2hiw': 'd', '1gwr': 'd', '2g24': 'd',
              '1c14': 'd', '1ywn': 'd', '1hvy': 'd', '1f9g': 'n', '1ai2': 'n', '2ivu': 'd',
              '2dq7': 'd', '1m2z': 'd', '2fb8': 'd', '1o5r': 'd', '2gh5': 'd', '1ke6': 'd',
              '1k7f': 'd', '1ucn': 'n', '1hw8': 'd', '2br1': 'd', '2i0e': 'd', '1js3': 'd',
              '1yqy': 'd', '1u4d': 'd', '1sqi': 'd', '2gsu': 'n', '1kvo': 'd', '1gpu': 'n',
              '1qpe': 'd', '1hvr': 'd', '1ig3': 'd', '1g7v': 'n', '1qmf': 'n', '1r58': 'd',
              '1v4s': 'd', '1fth': 'n', '1rsz': 'd', '1n2v': 'd', '1m17': 'd', '1kts': 'n',
              '1ywr': 'd', '2gyi': 'n', '1cg0': 'n', '5yas': 'n', '1icj': 'n', '1gkc': 'd',
              '1hqg': 'n', '1u30': 'd', '1nnc': 'n', '1c9y': 'n', '1j4i': 'd', '1qxo': 'n',
              '1o8b': 'n', '1nlj': 'n', '1rnt': 'n', '1d09': 'n', '1olq': 'n'}

    pdbs, tags = zip(*[[str(pdb), str(tag)] for pdb, tag in subset.items()])

    # with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
    #     dfs = executor.map(tractability_workflow, pdbs, tags)
    #
    # df = pd.concat(dfs, ignore_index=True)
    # df.sort_values(by='median', ascending=False)

    df = pd.read_csv('scores.csv')

    #
    df2 = pd.read_csv('/home/pcurran/covid/results/tractability/mpro.csv')
    df = pd.concat([df, df2])
    df.reset_index()
    #

    df = df.sort_values(by='median', ascending=False)
    joyplot(df, 'druggable_joy.png')

    # t = []
    # m = []
    # letter_to_number = {"d": 1, "n": 0}
    # for p in set(df['pdb']):
    #     a = df.loc[df['pdb'] == p]
    #     t.append(letter_to_number[list(a['tractability'])[0]])
    #     m.append(list(a['median'])[0])
    #
    # df = pd.DataFrame({'tractability': t,  'median':m})
    # df = df.sort_values(by='median', ascending=False)
    # data = list(zip(list(df["median"]), list(df["tractability"])))
    #
    # rocplot(data, fname='druggable_roc.png')


if __name__ == '__main__':
    main()
