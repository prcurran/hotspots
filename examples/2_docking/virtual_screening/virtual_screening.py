"""
Example script for the Hotspot API manuscript

The docking example has been adapted from the CCDC API documentation:
    - https://downloads.ccdc.cam.ac.uk/documentation/API/cookbook_examples/docking_examples.html
"""

import os
import gzip
import shutil
import tempfile
from shutil import copyfile
import concurrent

import numpy as np
import pandas as pd

from ccdc.docking import Docker
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.protein import Protein

from hotspots.calculation import Runner
from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.hs_docking import DockerSettings


def create_dataframe(fname):
    """
    Process the output ranking into a dataframe

    :param str fname: path to output file
    :return: `pandas.DataFrame`
    """
    lines = [l.strip("\n") for l in open(fname, "r").readlines()]

    header = lines[5]
    header = [b.strip() for b in [a for a in header.split("  ") if a != '' and a != '#']]

    print(header)

    all = lines[7:]
    cat = list(zip(*[[a for a in entry.split(" ") if a != ''] for entry in all]))
    df = pd.DataFrame({h: cat[i] for i, h in enumerate(header)})

    df["actives"] = np.array(map(lambda x: 'CHEMBL' in x, list(df['Ligand name'])))
    return df.sort_values(by=['Score'], ascending=False)


def dock(inputs):
    """


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
        copyfile(os.path.join(junk, "bestranking.lst"),
                 os.path.join(out_path, "bestranking.lst"))

    # GOLD docking routine
    ligand_path, out_path, hotspot, weight, search_efficiency = inputs
    docker = Docker()

    # GOLD settings
    docker.settings = DockerSettings()
    docker.settings.fitness_function = 'plp'
    docker.settings.autoscale = search_efficiency
    junk= os.path.join(out_path, "all")
    docker.settings.output_directory = junk

    # GOLD write lots of files we don't need in this example
    if not os.path.exists(junk):
        os.mkdir(junk)
    docker.settings.output_file = os.path.join(junk, "docked_ligands.mol2")

    # read the hotspot
    hotspot = HotspotReader(hotspot).read()
    for p, g in hotspot.super_grids.items():
        hotspot.super_grids[p] = g.max_value_of_neighbours()  # dilation to reduce noise

    add_ligands(docker,ligand_path)
    add_protein(docker, hotspot, junk)
    define_binding_site(docker, ligand_path)
    add_hotspot_constraint(docker, hotspot, weight)
    docker.dock(file_name=os.path.join(out_path, "hs_gold.conf"))
    write(docker, out_path)

    # Clean out unwanted files
    shutil.rmtree(junk)

    return f"docked {search_efficiency} {weight}"


def create_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)
    return path


def main():

    search_effiency = [10, 25, 50, 75, 100]
    weights = [10, 25, 50, 75, 100]

    parent = "/vagrant/github_pkgs/hotspots/examples/2_docking/virtual_screening/"

    ligand_path = os.path.join(parent, 'akt1')
    se = []
    wt = []
    otp = []

    results_path = create_dir(os.path.join(parent, "results"))
    for s in search_effiency:
        a = create_dir(os.path.join(results_path, f"search_{s}"))
        for w in weights:
            out_path = create_dir(os.path.join(a, f"weight_{w}"))

            otp.append(out_path)
            wt.append(w)
            se.append(s)

    ligand_path = [ligand_path] * len(se)
    hotspot = [os.path.join(parent, "hotspot", "out.zip")] * len(se)

    inputs = list(zip(ligand_path[:2], otp[:2], hotspot[:2], wt[:2], se[:2]))

    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        for z in executor.map(dock, inputs):
            print(z)

    for s in search_effiency:
        for w in weights:
            df = create_dataframe(os.path.join(results_path, f'search_{s}/weight_{w}/bestranking.lst'))
            df.to_csv(os.path.join(results_path, f'search_{s}/weight_{w}/bestranking.csv'))

if __name__ == "__main__":
    main()

