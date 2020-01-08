import concurrent
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from hotspots.calculation import Runner
from hotspots.result import Extractor


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
    for probe, grid in bcv_result.super_grids.items():
        values = grid.grid_values(threshold=5)
        median = np.median(values)

    # 4) return the data
    return pd.DataFrame({'scores': values,
                         'pdb': [protein] * len(values),
                         'median': [median] * len(values),
                         'tractability': [tag] * len(values),
                         })


def joyplot(df, fname='test.png'):
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

    palette = [(202 / 255, 41 / 255, 54 / 255),
               (140 / 255, 157 / 255, 196 / 255)]

    # Initialize the FacetGrid object
    ax = sns.FacetGrid(df, row="pdb", hue="tractability", height=4, aspect=75, size=.5, palette=palette)

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
    tag = {"d": "Druggable", "n": "Less-Druggable"}
    labels = [tag[s] for s in set(df["tractability"])]
    handles = [patches.Patch(color=col, label=lab) for col, lab in zip(palette, labels)]
    legend = plt.legend(handles=handles, title='Tractability', loc="upper right", bbox_to_anchor=(0.3, 7.5))
    frame = legend.get_frame()
    frame.set_facecolor('white')
    frame.set_edgecolor('white')

    plt.savefig(fname)
    plt.close()


def main():
    training_set = {'1e9x': 'd', '1hw8': 'd', '1sqi': 'd', '1r9o': 'd', '4cox': 'd', '1c14': 'd', '2bxr': 'd',
                    '2gh5': 'd', '1hvy': 'd', '1rsz': 'd', '1n2v': 'd', '1v4s': 'd', '1u4d': 'd', '1m17': 'd',
                    '2dq7': 'd', '1qpe': 'd', '1qhi': 'd', '2fb8': 'd', '1ke6': 'd', '2br1': 'd', '1ywr': 'd',
                    '2ivu': 'd', '2hiw': 'd', '2i0e': 'd', '1ywn': 'd', '1ig3': 'd', '1yvf': 'd', '1k8q': 'd',
                    '1kvo': 'd', '1xm6': 'd', '1udt': 'd', '1u30': 'd', '1r58': 'd', '1rwq': 'd', '1lpz': 'd',
                    '2g24': 'd', '1hvr': 'd', '1gkc': 'd', '1yqy': 'd', '1o5r': 'd', '1js3': 'd', '1k7f': 'd',
                    '1j4i': 'd', '1vbm': 'd', '1rv1': 'd', '1gwr': 'd', '1m2z': 'd', '3d4s': 'd', '1ai2': 'n',
                    '3pcm': 'n', '1d09': 'n', '1c9y': 'n', '1gpu': 'n', '1qmf': 'n', '1moq': 'n', '1ucn': 'n',
                    '1t03': 'n', '1qs4': 'n', '1fth': 'n', '1rnt': 'n', '1onz': 'n', '1x9d': 'n', '1nnc': 'n',
                    '1olq': 'n', '1jak': 'n', '1kts': 'n', '1nlj': 'n', '1icj': 'n', '1hqg': 'n', '2gsu': 'n',
                    '1g7v': 'n', '1f9g': 'n', '1qxo': 'n', '2gyi': 'n', '1o8b': 'n', '1cg0': 'n'}

    pdbs, tags = zip(*[[str(pdb), str(tag)] for pdb, tag in training_set.items()])

    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        dfs = executor.map(tractability_workflow, pdbs, tags)

    df = pd.concat(dfs, ignore_index=True)
    df = df.sort_values(by='median', ascending=False)

    df.to_csv('scores2.csv')
    joyplot(df, 'training_set.png')


if __name__ == '__main__':
    main()
