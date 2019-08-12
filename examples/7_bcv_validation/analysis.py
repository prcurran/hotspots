"""

generate paper figures

"""
from __future__ import print_function

import ast
import datetime
import os

import numpy as np
import pandas as pd
from ccdc.cavity import Cavity
from ccdc.io import MoleculeReader
from pipeline import HotspotPipeline

from hotspots.hs_io import HotspotReader
from hotspots.grid_extension import Grid


class Hot(HotspotPipeline):
    def _get_ligand_cavity(self):
        self.ligand_cavity = os.path.join(self.working_dir, "ligand_cavity.dat")
        tolerance = 0
        lc = []
        mols = [MoleculeReader(path)[0] for other_id, lig_dic in self.extracted_ligands.items()
                for l, path in lig_dic.items()]

        point = [round(np.mean([a.coordinates.x for mol in mols for a in mol.heavy_atoms])),
                 round(np.mean([a.coordinates.y for mol in mols for a in mol.heavy_atoms])),
                 round(np.mean([a.coordinates.z for mol in mols for a in mol.heavy_atoms]))]

        cavs = Cavity.from_pdb_file(self.apo_prep)

        for i, c in enumerate(cavs):

            mini = c.bounding_box[0]
            maxi = c.bounding_box[1]
            if all([mini.x - tolerance < point[0] < maxi.x + tolerance,
                    mini.y - tolerance < point[1] < maxi.y + tolerance,
                    mini.z - tolerance < point[2] < maxi.z + tolerance]):

                lc.append(int(i))

        with open(self.ligand_cavity, "w") as f:
            for cav_id in lc:
                f.write(str(cav_id))

        if len(lc) < 1:
            print("no ligand cavity for {}".format(self.apo))
            return lc

        else:
            return lc

    def _get_overlap_dic(self, cav, other_id, lig_id):
        """


        :param cav:
        :param other_id:
        :param lig_id:
        :return:
        """
        if os.path.exists(self.atomic_overlaps[cav][other_id][lig_id]):
            with open(self.atomic_overlaps[cav][other_id][lig_id], "r") as f:
                x = f.read()
            return ast.literal_eval(x)
        else:
            return None

    def _get_vol_data(self, cav, other_id, lig_id):
        """
        extract the volume data

        :param int cav: cavity identifier
        :param str other_id: protein identifier
        :param str lig_id: ligand identifier
        :return:
        """
        if os.path.exists(self.all_overlaps[cav][other_id][lig_id]):
            with open(self.all_overlaps[cav][other_id][lig_id], "r") as f:
                vol = f.read()
            return float(vol)
        else:
            return 0.0

    def _is_top_cavity(self, cav):
        """
        determine if top ranked cavity

        :param cav:
        :return:
        """
        if os.path.exists(self.cavity_rank):
            with open(self.cavity_rank, "r") as f:
                print(self.cavity_rank)
                c = f.read()
                try:
                    if int(cav) == int(c):
                        return True
                    else:
                        return False
                except:
                    return "NaN"
        else:
            return "NaN"

    def _max_overlap(self, l):
        try:
            vals = []
            for cav_id, prot_dic in self.all_overlaps.items():
                for prot_id, lig_dic in prot_dic.items():
                    for lig_id, path in lig_dic.items():
                        if os.path.exists(path) and l == lig_id:
                            with open(path, 'r') as f:
                                vals.append(float(f.read()))
            return max(vals)
        except:
            return 0

    def _get_cavity_score(self, cav_id):
        """
        retrieve cavity score

        :param cav_id:
        :return:
        """
        if os.path.exists(self.cavity_score[cav_id]):
            with open(self.cavity_score[cav_id], 'r') as f:
                score = float(f.read())
            return score
        else:
            return "NaN"

    def analysis(self):
        """
        report the volume analysis

        :return:
        """
        keys = ["apo", "buriedness_method", "cavity_id", "other_id", "ligand_id", "cavity_score",
                "top_cavity", "volume_overlap", "atomic_label", "atomic_overlap", "atom_type", "ligand_cavity"]
        data = []
        self._get_cavities(min_vol=200)
        lc = self._get_ligand_cavity()

        for cav in range(len(self.cavities)):
            v1 = self._is_top_cavity(cav)
            score = self._get_cavity_score(cav)
            if int(cav) in lc:
                lc_bool = True
            else:
                lc_bool = False
            for i, prot_id in enumerate(self.protein_id):
                for lig_id in self.ligand_id[i]:
                    atomic_overlap_dic = self._get_overlap_dic(cav, prot_id, lig_id)
                    v2 = self._get_vol_data(cav, prot_id, lig_id)

                    if atomic_overlap_dic == None:
                        pass
                    else:
                        for n, atom_dic in atomic_overlap_dic.items():
                            for label, overlap in atom_dic.items():
                                data.append(
                                    [self.apo, self.buriedness_method, cav, prot_id, lig_id, score, v1, v2,
                                     label, overlap, n, lc_bool])

        return pd.DataFrame(dict(zip(keys, (zip(*data)))))

    def table(self):
        """
        report the volume analysis

        :return:
        """
        keys = ["apo", "prot_id", "ligand_id", self.buriedness_method]
        data = []
        self._get_cavities(min_vol=200)
        lc = self._get_ligand_cavity()
        if len(lc) == 0:
            for i, prot_id in enumerate(self.protein_id):
                for lig_id in self.ligand_id[i]:
                    v2 = self._get_vol_data(0, prot_id, lig_id)
                    data.append([self.apo, prot_id, lig_id, v2])

        else:
            for cav in range(len(self.cavities)):
                if int(cav) in lc:
                    for i, prot_id in enumerate(self.protein_id):
                        for lig_id in self.ligand_id[i]:
                            v2 = self._get_vol_data(cav, prot_id, lig_id)
                            data.append([self.apo, prot_id, lig_id, v2])

        return pd.DataFrame(dict(zip(keys, (zip(*data)))))

    def bcv_effect(self):
        pdb = []
        hot = []
        bcv = []
        self._get_cavities(min_vol=200)
        for cav in range(len(self.cavities)):
            lc = self._get_ligand_cavity()
            print(lc)
            if int(cav) in lc:
                print(self.hotspot[cav])
                ahr = HotspotReader(os.path.join(self.hotspot[cav], "out.zip")).read()
                a = Grid.get_single_grid(ahr.super_grids, mask=False)
                for i, prot_id in enumerate(self.protein_id):
                    for lig_id in self.ligand_id[i]:
                        bhr = HotspotReader(os.path.join(self.bcv[cav][prot_id][lig_id], "out.zip")).read()
                        b = Grid.get_single_grid(bhr.super_grids, mask=False)
                        pdb.append(self.apo)
                        hot.append(a.count_grid())
                        bcv.append(b.count_grid())

        return pdb, hot, bcv


def main():
    prefix = "/vagrant/github_pkgs/hotspots/examples/7_bcv_validation"
    buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']
    df = pd.read_csv("inputs.csv")
    frags = set(df['fragment'])
    leads = set(df['lead'])

    hot_pdbs = set(df['apo'])
    reports = []

    for i, pdb in enumerate(hot_pdbs):
        for method in buriedness_methods:
            ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
            proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])

            hp = Hot(apo=pdb, buriedness_method=method, protein_id=proteins, ligand_id=ligands)
            report = hp.analysis()
            reports.append(report)
            print(report)

    dat = pd.concat(reports, ignore_index=True)

    classification = []
    for a in list(dat['other_id']):
        if a in frags:
            classification.append("fragment")
        elif a in leads:
            classification.append("lead")

    dat["lig_class"] = classification
    dat.to_csv("analysis.csv")

def table():
    buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']
    df = pd.read_csv("inputs.csv")
    hot_pdbs = set(df['apo'])
    frags = set(df['fragment'])
    leads = set(df['lead'])
    reports = []

    for i, pdb in enumerate(hot_pdbs):
        to_merge = []
        for method in buriedness_methods:
            ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
            proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])

            hp = Hot(apo=pdb, buriedness_method=method, protein_id=proteins, ligand_id=ligands)
            report = hp.table()
            to_merge.append(report)
        z = pd.merge(to_merge[0], to_merge[1], how='outer', on=['apo', 'prot_id', 'ligand_id'])
        reports.append(pd.merge(z, to_merge[2], how='outer', on=['apo', 'prot_id', 'ligand_id']))


    dat = pd.concat(reports, ignore_index=True)
    
    frags = set(df['fragment'])
    leads = set(df['lead'])
    classification = []
    for a in list(dat['other_id']):
        if a in frags:
            classification.append("fragment")
        elif a in leads:
            classification.append("lead")

    dat["lig_class"] = classification
    dat.to_csv("table.csv")


def bcv_effect():
    print("here")
    method = 'ghecom'
    df = pd.read_csv("inputs.csv")
    hot_pdbs = set(df['apo'])
    print(hot_pdbs)
    p = []
    h = []
    b = []
    for i, pdb in enumerate(list(hot_pdbs)):
        ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
        proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])

        hp = Hot(apo=pdb, buriedness_method=method, protein_id=proteins, ligand_id=ligands)
        print(hp)
        pdbcode, hot, bcv = hp.bcv_effect()
        p.extend(pdbcode)
        h.extend(hot)
        b.extend(bcv)

    f = pd.DataFrame({"pdb": p, "hotspot": h, "bcv": b})
    f.to_csv("bcv_effect.csv")

if __name__ == "__main__":
    # table()
    bcv_effect()
