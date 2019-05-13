"""

This script wraps around the Protein Plus web services
    - SIENA: automated construction of protein ensembles from the PDB

    - Bietz, S. Rarey, M.: SIENA: Efficient Compilation of Selective Protein Binding Site Ensembles.
    Journal of Chemical Information and Modeling,56(1): 248-59.

    - Bietz, S. Rarey, M.: ASCONA: Rapid Detection and Alignment of Protein Binding Site Conformations.
    Journal of Chemical Information and Modeling, 55(8):1747-1756.

"""
import shlex, subprocess, ast, time, json, os, csv, re, urllib2, urlparse, tempfile
from collections import OrderedDict
import numba


@numba.njit()
def tanimoto_dist(a, b):
    """
    calculate the tanimoto distance between two fingerprint arrays
    :param a:
    :param b:
    :return:
    """
    dotprod = np.dot(a, b)
    tc = dotprod / (np.sum(a) + np.sum(b) - dotprod)
    return 1.0 - tc


class Ensemble(object):
    """
    TO DO: process the output data here (create ligands, proteins properties etc)

    Object containing data on the Ensemble members for a particular target

    Outputted data:

        - result_table - main result, binding site ensemble with different conformations of the same/closely
        related binding sites found in the PDB (CSV-table, one PDB-entry per row)
        - pdb_files - alternative conformations of binding site from the PDB (list of PDB-files)
        - ligands - ligands at the binding site of the alternative conformations (list of SDF-files,
        positions match the PDB entries in the pdb_files-list)
        - alignment - shows positions in the different chains of the found PDB-entries that can be aligned (TXT-file)
        - parameters - the values of parameters used for this request

    """
    def __init__(self, data):
        self.data = data
        self.to_retrieve = ['result_table', 'pdb_files', 'ligands', 'alignment']
        self._ligands = None


    def _write(self, url, out_dir):
        """
        handles writing all the data

        :param str url: url for aligned files
        :param str out_dir: path of output directory
        """
        print url
        url = re.sub('/esults/', '/results/', url)              # fix url
        a = urlparse.urlparse(url)
        fname = os.path.basename(a.path)
        ext = fname.split(".")[1]
        path = os.path.join(out_dir, fname)

        with open(path, "wb") as w:
            if ext == "csv":
                # format csv
                stri = [[a for a in line.split(";") if a != '\n'] for line in urllib2.urlopen(url).readlines()]
                csv_writer = csv.writer(w, delimiter=',')
                for s in stri:
                    csv_writer.writerow(s)
            else:
                # otherwise no formatting requried
                w.write(urllib2.urlopen(url).read())

    def save(self, out_dir):
        """
        save the ensemble data to output directory

        :param str out_dir: path to output directory
        """
        for t in self.to_retrieve:
            new_out_dir = os.path.join(out_dir, t)
            if not os.path.exists(new_out_dir):
                os.mkdir(new_out_dir)
            urls = self.data[t]
            if type(urls) is list:
                for url in urls:
                    self._write(url, new_out_dir)
            else:
                self._write(urls, new_out_dir)

    # @staticmethod
    # def _cluster_ligands(ligands, t):
    #     """
    #
    #     :return:
    #     """
    #     def fingerprint_array(ligands):
    #         X =[]
    #         for l in ligands:
    #             arr = np.zeros((0,))
    #             fp = AllChem.GetMorganFingerprintAsBitVect(l, 2)
    #             DataStructs.ConvertToNumpyArray(fp, arr)
    #             X.append(arr)
    #         return X
    #
    #     cluster_dic = {}
    #
    #     # generate fingerprint array
    #     X = fingerprint_array(ligands)
    #     if len(X) < 2:
    #         X = fingerprint_array(ligands)
    #         if len(X) < 2:
    #             raise ValueError("Fingerprint array must contain more than 1 entry")
    #
    #     # dimensionality reduction
    #     tsne_X = TSNE(n_components=2, metric=tanimoto_dist).fit_transform(np.array(X, dtype=np.float32))
    #
    #     # clustering
    #     cluster_tsne = hdbscan.HDBSCAN(min_cluster_size=2, gen_min_span_tree=True)
    #     cluster_tsne.fit(tsne_X)
    #
    #     for i, label in enumerate(cluster_tsne.labels_):
    #         if label == -1:
    #             continue
    #         else:
    #             if label in cluster_dic:
    #                 cluster_dic[label].append(ligands[i])
    #             else:
    #                 cluster_dic.update({label: [ligands[i]]})
    #
    #     x = [tsne_X.T[0][j] for j, l in enumerate(cluster_tsne.labels_) if l != -1]
    #     y = [tsne_X.T[1][j] for j, l in enumerate(cluster_tsne.labels_) if l != -1]
    #     hue = [l for j, l in enumerate(cluster_tsne.labels_) if l != -1]
    #
    #     plt.scatter(x, y, c=hue, cmap='RdBu')
    #
    #     plt.title("{} clusters".format(t))
    #     plt.savefig("{}.png".format(t))
    #     plt.close()
    #     if len(cluster_dic) == 0:
    #         print("NO CLUSTERS FOUND")
    #         cluster_dic = {i: [ligands[i]] for i in range(0, len(ligands))}
    #
    #     return cluster_dic

    # def select_diverse_ligands(self, target):
    #     """
    #
    #
    #     :return:
    #     """
    #     # download ligands
    #     # convert to fps
    #     # cluster on ligands
    #     # update data to only include selected members
    #     tmp = tempfile.mkdtemp()
    #     ligands = self.data['ligands']
    #     for ligand in ligands:
    #         self._write(ligand, tmp)
    #
    #     files = [os.path.join(tmp, f) for f in os.listdir(tmp) if os.path.isfile(os.path.join(tmp, f))]
    #     ligands = {os.path.basename(f).split(".")[0]: x
    #                for f in files for x in Chem.ForwardSDMolSupplier(f) if x is not None}
    #
    #     for n, l in ligands.items():
    #         l.SetProp("_Name", n)
    #
    #     cluster_dict = self._cluster_ligands(ligands=ligands, t=target)
    #     reps = [l[0] for l in cluster_dict.values() if len(l) != 0]
    #
    #     print reps


class Search(object):
    """
    post request to initiate job
    get request to retrieve data (iteratively until the job has complete)

    >>> searcher = Search()
    >>> ensemble = searcher.create_ensemble(pdb_code="1aq1",
                                        ligand="STU_A_299",
                                        reduction_procedure="",
                                        num_members="")

    >>> ensemble.save(out_dir = "/home/pcurran/patel/CDK2/1aq1_ensemble")

    """
    class Settings(object):
        """


        """
        def __init__(self):
            self.url = 'https://proteins.plus/api/siena_rest'
            self.data = {"siena": {
                            "pdbCode":"",
                            "pocket":"",
                            "ligand":"",
                            "mode":""}
                         }
            # self.data = {"reduction_procedure":"",
            #              "siena": {"bb_clustering":"",
            #                        "all_atom_clustering":"",
            #                        "ligand_driven_selection":"",
            #                        "ligand": "",
            #                        "pocket":"",
            #                        "pdbCode": "",
            #                        "siteRadius":"6.5",
            #                        "fragment_length": "10",
            #                        "flexibility_sensitivity": "0.6",
            #                        "fragment_distance": "4",
            #                        "minimalSiteIdentity": "0.7",
            #                        "minimalSiteCoverage":"",
            #                        "maximum_mutations":"",
            #                        "resolution":"",
            #                        "maximumBackbone":"",
            #                        "depositionYear":"",
            #                        "ecNumber":"",
            #                        "electronDensityAvailable": "",
            #                        "identicalGlobalSequence": "",
            #                        "noGlobalMutations": "true",
            #                        "unique_sequence":"",
            #                        "holo_only": "true",
            #                        "unique_ligands": "",
            #                        "complete_residues_only": ""
            #                        }
            #              }

    def __init__(self, settings=None):
        """
        Search initialisation

        :param `Search.Settings` settings: a settings object, modify to change advance settings
        """
        if settings is None:
            self.settings = self.Settings()
        else:
            self.settings = settings

    def create_ensemble(self, pdb_code, ligand, mode, pocket=None):
        """


        :param str pdb_code: PDB code to base the ensemble on
        :param str ligand: ligand identifier (example format: STU_A_299)
        :param str reduction_procedure: either bb_clustering, all_atom_clustering, ligand_driven_selection
        :param str num_members: number of members in the ensemble
        :return:
        """
        self.settings.data['siena']["pdbCode"] = pdb_code
        if pocket is not None:
            self.settings.data['siena']['pocket'] = pocket
        self.settings.data['siena']["ligand"] = ligand
        self.settings.data['siena']['mode'] = mode

        results_url = self._post()
        return Ensemble(self._get(results_url))

    def _run(self, cmd):
        """
        runs commandline procedure (urllib doesn't work for some reason)

        :param str cmd: command line str using curl
        :return:
        """
        args = shlex.split(cmd)
        return subprocess.check_output(args)

    def _post(self):
        """
        Initiate the job using a POST request

        :return:
        """
        cmd = """curl -d '{}' -H "Accept: application/json" -H "Content-Type: application/json" -X POST {}"""\
            .format(json.dumps(self.settings.data, ensure_ascii=True, indent=2, default=True),self.settings.url)

        response = ast.literal_eval(self._run(cmd))
        print response
        return response['location']

    def _get(self, results_url):
        """
        Collect the results using GET. May have to try multiple times.

        :return:
        """
        cmd = """curl {}""".format(results_url)

        response = ast.literal_eval(self._run(cmd))
        if response["status_code"] == 400:
            print "status code: {}".format(response['status_code'])
            raise RuntimeError()

        status_code = int(response["status_code"])
        while status_code > 200:
            print "status code: {}".format(response['status_code'])
            time.sleep(10)
            response = ast.literal_eval(self._run(cmd))
            status_code = int(response["status_code"])
            if status_code == 400:
                raise RuntimeError("help")

        # create Ensemble

        print type(response["status_code"])
        return response


def main():
    searcher = Search()
    ensemble = searcher.create_ensemble(pdb_code="1aq1",
                                        ligand="STU_A_299",
                                        mode="ligand_pose_comparison")

    ensemble.save(out_dir = "/home/pcurran/patel/CDK2/1aq1_ensemble")
    #ensemble.select_diverse_ligands(target="CDK2")

if __name__ == "__main__":
    main()

