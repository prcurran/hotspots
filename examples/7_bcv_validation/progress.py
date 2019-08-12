"""

track the progress of the BCV validation

"""
from __future__ import print_function
import datetime
import pandas as pd
from pipeline import HotspotPipeline


class Hot(pipeline.HotspotsPipeline):
    def _cav_volume(self, cav):
        if os.path.exists(self.cavities_volumes[cav]):
            with open(self.cavities_volumes[cav], "r") as f:
                vol = f.read()
            return float(vol)
        else:
            return "NaN"

    def _step_1_status(self):
        if os.path.exists(self.apo_prep):
            return 1
        else:
            return 0

    def _step_2_status(self, other_id, lig_id):
        if os.path.exists(self.extracted_ligands[other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_3_status(self, other_id, lig_id):
        if os.path.exists(self.ligand_volume[other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_4_status(self):
        if os.path.exists(self.buriedness):
            return 1
        else:
            return 0

    def _step_5_status(self, cav):
        if os.path.exists(self.cavities[cav]):
            return 1
        else:
            return 0

    def _step_6_status(self, cav):
        if os.path.exists(os.path.join(self.superstar[cav], "out.zip")):
            return 1
        else:
            return 0

    def _step_7_status(self, cav):
        if os.path.exists(os.path.join(self.hotspot[cav], "out.zip")):
            return 1
        else:
            return 0

    def _step_8_status(self, cav, other_id, lig_id):
        if os.path.exists(os.path.join(self.bcv[cav][other_id][lig_id], "out.zip")):
            return 1
        else:
            return 0

    def _step_9_status(self, cav, other_id, lig_id):
        if os.path.exists(self.all_overlaps[cav][other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_10_status(self, cav, other_id, lig_id):
        if os.path.exists(self.matched[cav][other_id][lig_id]):
            return 1
        else:
            return 0

    def _step_11_status(self):
        if os.path.exists(self.cavity_rank):
            return 1
        else:
            return 0

    def status(self):
        """
        report the status of the calculations

        :return:
        """
        keys = ["apo", "buriedness_method", "cavity_id", "other_id", "ligand_id", "cav_vol", "apo_prep", "extract_ligand",
                "extract_vol", "calculate_buriedness", "cavity_detection", "calculate_superstar", "calculate_hotspots", "calculate_bcv",
                "volume_overlap", "matched_atoms", "cavity_rank"]
        data = []
        self._get_cavities(min_vol=200)

        s1 = self._step_1_status()
        s4 = self._step_4_status()
        s11 = self._step_11_status()

        for cav in range(len(self.cavities)):
            vol = self._cav_volume(cav)
            s5 = self._step_5_status(cav)
            s6 = self._step_6_status(cav)
            s7 = self._step_7_status(cav)
            for i, prot_id in enumerate(self.protein_id):
                for lig_id in self.ligand_id[i]:
                    s2 = self._step_2_status(prot_id, lig_id)
                    s3 = self._step_3_status(prot_id, lig_id)
                    s8 = self._step_8_status(cav, prot_id, lig_id)
                    s9 = self._step_9_status(cav, prot_id, lig_id)
                    s10 = self._step_10_status(cav, prot_id, lig_id)

                    data.append([self.apo, self.buriedness_method, cav, prot_id, lig_id, vol, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11])

        return pd.DataFrame(dict(zip(keys, (zip(*data)))))


def main():
    prefix = "/vagrant/github_pkgs/hotspots/examples/7_bcv_validation"
    buriedness_methods = ['ligsite', 'ghecom', 'ghecom_internal']
    df = pd.read_csv("inputs.csv")
    hot_pdbs = set(df['apo'])
    reports = []

    for i, pdb in enumerate(hot_pdbs):
        for method in buriedness_methods:
            ligands = list(df.loc[df['apo'] == pdb]['fragment_ID']) + list(df.loc[df['apo'] == pdb]['lead_ID'])
            proteins = list(df.loc[df['apo'] == pdb]['fragment']) + list(df.loc[df['apo'] == pdb]['lead'])

            hp = Hot(apo=pdb, buriedness_method=method, protein_id=proteins, ligand_id=ligands)
            report = hp.status()
            reports.append(report)

    progress = pd.concat(reports)

    date = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    progress.to_csv("progress_report_{}.csv".format(date))
    print(progress)


if __name__ == "__main__":
    main()