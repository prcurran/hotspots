"""
Example script for the Hotspot API manuscript

The docking example has been adapted from the CCDC API documentation:
    - https://downloads.ccdc.cam.ac.uk/documentation/API/cookbook_examples/docking_examples.html
"""

import gzip
import os
import shutil
from shutil import copyfile

import numpy as np
import operator
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ccdc.descriptors import StatisticalDescriptors as sd
from ccdc.docking import Docker
from ccdc.io import MoleculeReader, MoleculeWriter
from hotspots.hs_docking import DockerSettings
from hotspots.hs_io import HotspotReader



def dock(inputs):
    """
    submit a GOLD API docking calculation using docking constraints automatically generated from the Hotspot API

    :param ligand_path:
    :param out_path:
    :param hotspot:
    :param weight:
    :return:
    """

    def add_ligands(docker, ligand_path):

        with gzip.open(os.path.join(ligand_path, "actives_final.mol2.gz"), 'rb') as f_in:
            with open(os.path.join(docker.settings.output_directory, "actives_final.mol2"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        with gzip.open(os.path.join(ligand_path, "decoys_final.mol2.gz"), 'rb') as f_in:
            with open(os.path.join(docker.settings.output_directory, "decoys_final.mol2"), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        docker.settings.add_ligand_file(os.path.join(docker.settings.output_directory,
                                                     "actives_final.mol2"),
                                        ndocks=5)

        docker.settings.add_ligand_file(os.path.join(docker.settings.output_directory,
                                                     "decoys_final.mol2"),
                                        ndocks=5)

    def add_protein(docker, hotspot, junk):

        pfile = os.path.join(junk, "protein.mol2")
        with MoleculeWriter(pfile) as w:
            w.write(hotspot.protein)

        docker.settings.add_protein_file(pfile)

    def define_binding_site(docker, ligand_path):

        crystal_ligand = MoleculeReader(os.path.join(ligand_path, "crystal_ligand.mol2"))[0]
        docker.settings.binding_site = docker.settings.BindingSiteFromLigand(protein=docker.settings.proteins[0],
                                                                             ligand=crystal_ligand)

    def add_hotspot_constraint(docker, hotspot, weight):

        if int(weight) != 0:
            constraints = docker.settings.HotspotHBondConstraint.create(protein=docker.settings.proteins[0],
                                                                        hr=hotspot,
                                                                        weight=int(weight),
                                                                        min_hbond_score=0.05,
                                                                        max_constraints=1)

            for constraint in constraints:
                docker.settings.add_constraint(constraint)

    def write(docker, out_path):

        results = Docker.Results(docker.settings)

        # write ligands
        with MoleculeWriter(os.path.join(out_path, "docked_ligand.mol2")) as w:
            for d in results.ligands:
                w.write(d.molecule)

        # copy ranking file
        # in this example, this is the only file we use for analysis. However, other output files can be useful.
        copyfile(os.path.join(junk, "bestranking.lst"),
                 os.path.join(out_path, "bestranking.lst"))

    # GOLD docking routine
    ligand_path, out_path, hotspot, weight, search_efficiency = inputs
    docker = Docker()

    # GOLD settings
    docker.settings = DockerSettings()
    docker.settings.fitness_function = 'plp'
    docker.settings.autoscale = search_efficiency
    junk = os.path.join(out_path, "all")
    docker.settings.output_directory = junk

    # GOLD write lots of files we don't need in this example
    if not os.path.exists(junk):
        os.mkdir(junk)
    docker.settings.output_file = os.path.join(junk, "docked_ligands.mol2")

    # read the hotspot
    hotspot = HotspotReader(hotspot).read()
    for p, g in hotspot.super_grids.items():
        hotspot.super_grids[p] = g.max_value_of_neighbours()  # dilation to reduce noise

    add_ligands(docker, ligand_path)
    add_protein(docker, hotspot, junk)
    define_binding_site(docker, ligand_path)
    add_hotspot_constraint(docker, hotspot, weight)
    docker.dock(file_name=os.path.join(out_path, "hs_gold.conf"))
    write(docker, out_path)

    # Clean out unwanted files
    shutil.rmtree(junk)


def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
    return path


def rank_stats(parent, s, w):
    """
    Create two `pandas.DataFrame` from the "bestranking.lst"
    GOLD output

    :param str parent: path to parent directory
    :param str s: name of subdirectory
    :param str w: name of subsubdirectory
    :return:
    """
    # read data and process data from output file
    fname = os.path.join(parent,
                         f"search_efficiency_{s}",
                         str(w),
                         "bestranking.lst")
    lines = [l.strip("\n") for l in open(fname, "r").readlines()]
    header = lines[5]
    header = [b.strip() for b in
              [a for a in header.split("  ") if a != '' and a != '#']]
    all = lines[7:]

    cat = list(zip(*[[a for a in entry.split(" ") if a != '']
                     for entry in all]))

    # generate a dataframe and alter datatypes
    df = pd.DataFrame({h: cat[i] for i, h in enumerate(header)})
    df["actives"] = np.array(
        list(map(lambda x: 'CHEMBL' in x, list(df['Ligand name'])))
    ).astype(int)
    df["search_efficiency"] = [int(s)] * len(df)
    df["weight_int"] = [int(w)] * len(df)
    df["weight_str"] = [str(w)] * len(df)
    df["Score"] = df["Score"].astype(float)
    df["time"] = df["time"].astype(float)
    df["log_time"] = np.log10(list(df["time"]))
    df = df[['Score',
             'log_time',
             'actives',
             'search_efficiency',
             'weight_int',
             'weight_str']]
    df = df.sort_values(by=['Score'], ascending=False)

    # Use CCDC's descriptors API
    rs = sd.RankStatistics(scores=list(zip(list(df['Score']),
                                           list(df['actives']))),
                           activity_column=operator.itemgetter(1))

    # ROC
    tpr, fpr = rs.ROC()
    df["tpr"] = tpr
    df["fpr"] = fpr

    # Enrichment Metrics
    metric_df = pd.DataFrame({"search efficiency": [s],
                             "weight": [w],
                             "AUC": [rs.AUC()],
                             "EF1": [rs.EF(fraction=0.01)],
                             "EF5": [rs.EF(fraction=0.05)],
                             "EF10": [rs.EF(fraction=0.1)],
                             "BEDROC16": [rs.BEDROC(alpha=16.1)],
                             "BEDROC8": [rs.BEDROC(alpha=8)]
                             })
    return df, metric_df


def roc_plot(df, search_efficiency, ax, palette='Set2'):
    """
    Plot a ROC plot for the docking data

    :param `pands.DataFrame` df: data
    :param int search_efficiency: data is split by search efficiency
    :param `matplotlib.axes.Axes` ax: Matplotlib axis to plot data onto
    :param list palette: list of RGB tuples
    :return:
    """

    selected = df.loc[df['search_efficiency'] == search_efficiency]
    # random
    sns.lineplot(x=[0, 1], y=[0, 1], color="grey", ax=ax)
    # docking rank

    d = {"0": {"color": sns.color_palette(palette)[0],"linestyle": "-"},
         "10": {"color": sns.color_palette(palette)[1],  "linestyle": "-"},
         "100": {"color": sns.color_palette(palette)[2], "linestyle": "-"}}

    lines = [ax.plot(grp.fpr, grp.tpr, label=n, **d[n])[0]
             for n, grp in selected.groupby("weight_str")]


    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticks([0, 0.5, 1])
    ax.set_title(label=f"Search Efficiency = {search_efficiency}%",
                 fontdict={'fontsize':10})

    return lines


def box_plot(df, search_efficiency, ax, palette='Set2'):
    """
    Plot a boxplot for the docking time data

    :param `pands.DataFrame` df: data
    :param int search_efficiency: data is split by search efficiency
    :param `matplotlib.axes.Axes` ax: Matplotlib axis to plot data onto
    :param list palette: list of RGB tuples
    :return:
    """

    selected = df.loc[df['search_efficiency'] == search_efficiency]
    sns.boxplot(x="weight_int", y="log_time", order=[0,10,100],data=selected,
                palette=palette,linewidth=1.2, fliersize=0.5, ax=ax)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticks([])

    # just for asthetics
    for patch in ax.artists:
        r,g,b,a = patch.get_edgecolor()
        patch.set_edgecolor((r,g,b,.1))


def _asthetics(fig, axs, lines):
    """
    Extra formatting tasks

    :param `matplotlib.figure.Figure` fig: mpl figure
    :param list axs: list of  `matplotlib.axes.Axes`
    :param `matplotlib.lines.Line2D` lines: ROC lines
    :return:
    """
    # format axes
    axs[0][0].set_yticks([0, 0.5, 1])
    yticks = [-1, 0, 1, 2]
    axs[1][0].set_yticks(yticks)
    axs[1][0].set_yticklabels([str(10 ** float(l)) for l in yticks])
    axs[0][1].set_xlabel("FPR")
    axs[0][0].set_ylabel("TPR")
    axs[1][0].set_ylabel("Time per Molecule (s)")
    # format legend
    fig.legend(lines,
               [f"{w}" for w in [0,10,100]],
               (.83, .42),
               title="Constraint Weight")
    # format canvas
    plt.subplots_adjust(left=0.1, right=0.8, top=0.86)

def main():
    # read and format the data
    search_effiencies = [1, 10, 100]
    weights = [0, 10, 100]
    parent = "/home/pcurran/github_packages/hotspots/examples/2_docking/virtual_screening/akt1"

    ligand_path = parent
    out_path = os.path.join(parent, "everything")
    hotspot = os.path.join(parent, "hotspot", "out.zip")
    weight = 0
    search_efficiency = 1

    inputs = (ligand_path, out_path, hotspot, weight, search_efficiency)

    dock(inputs)

    # df1, df2 = zip(*[rank_stats(parent, s, w)
    #                  for s in search_effiencies for w in weights])
    #
    # # Plotted Data (ROC and Box plots)
    # df1 = pd.concat(df1, ignore_index=True)
    #
    # # Table Data (Rank Stats: AUC, EF, BEDROC)
    # df2 = pd.concat(df2, ignore_index=True)
    # df2.to_csv('rankstats.csv')
    #
    # # Plot the ROC and box plots
    # sns.set_style('white')
    # fig, axs = plt.subplots(nrows=2,
    #                         ncols=3,
    #                         sharey='row',
    #                         gridspec_kw={'wspace':0.26,
    #                                      'hspace':0.22},
    #                         figsize=(10,6), dpi=200)
    #
    # for i, row in enumerate(axs):
    #     for j, ax in enumerate(row):
    #         if i == 0:
    #             lines = roc_plot(df1, search_effiencies[j], ax)
    #         else:
    #             box_plot(df1, search_effiencies[j], ax)
    #
    # _asthetics(fig, axs, lines)
    # plt.savefig("new_grid.png")
    # plt.close()

if __name__ == "__main__":
    main()
