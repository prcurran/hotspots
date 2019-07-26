# import shlex, subprocess, ast, time, json, os, csv, re, urllib, urlparse, tempfile
# from ccdc.protein import Protein
# from ccdc.io import MoleculeReader
#
#
# class Result(object):
#
#     def __init__(self, data, out_dir):
#         self.data = data
#         self.files = {'protein':[], 'ligands':[]}
#         self.out_dir = out_dir
#         self._write()
#         self._protein = Protein.from_file(self.files['protein'][0])
#         self._ligands = [y for y in [MoleculeReader(x) for x in self.files['ligands']]]
#
#
#     def _write(self):
#         """
#         save the ensemble data to output directory
#
#         :param str out_dir: path to output directory
#         """
#         for t in self.files.keys():
#             new_out_dir = os.path.join(self.out_dir, t)
#             if not os.path.exists(new_out_dir):
#                 os.mkdir(new_out_dir)
#
#             urls = self.data[t]
#
#             if type(urls) is list:
#                 for url in urls:
#                     self.files[t].append(self._out(url, new_out_dir))
#             else:
#                 self.files[t].append(self._out(urls, new_out_dir))
#
#     def _out(self, url, out_dir):
#         """
#
#         :param url:
#         :param out_dir:
#         :return:
#         """
#         url = re.sub('/esults/', '/results/', url)              # fix url
#         a = urlparse.urlparse(url)
#         fname = os.path.basename(a.path)
#         str = urllib2.urlopen(url).read()
#         path = os.path.join(out_dir, fname)
#         with open(path, "wb") as w:
#             w.write(str)
#         return path
#
#     @property
#     def protein(self):
#         return self._protein
#
#     @property
#     def ligands(self):
#         return self._ligands
#
#
# class Protoss(object):
#     """
#
#     """
#     def __init__(self, out_dir=tempfile.mkdtemp()):
#         """
#         """
#         self.url = "https://proteins.plus/api/protoss_rest"
#         self.data = {"protoss": {"pdbCode": ""}}
#         self.pdb = ""
#         self.out_dir = out_dir
#
#     def add_hydrogens(self, pdb_code):
#         """
#
#         """
#         self.pdb = pdb_code
#         self.data['protoss']["pdbCode"] = pdb_code
#
#         results_url = self._post()
#         return Result(data=self._get(results_url), out_dir=self.out_dir)
#
#     def _run(self, cmd):
#         """
#         runs commandline procedure (urllib doesn't work for some reason)
#
#         :param str cmd: command line str using curl
#         :return:
#         """
#         args = shlex.split(cmd)
#         return subprocess.check_output(args)
#
#     def _post(self):
#         """
#         Initiate the job using a POST request
#
#         :return:
#         """
#         cmd = """curl -d '{}' -H "Accept: application/json" -H "Content-Type: application/json" -X POST {}"""\
#             .format(json.dumps(self.data, ensure_ascii=True, indent=2, default=True),self.url)
#
#         response = ast.literal_eval(self._run(cmd))
#         return response['location']
#
#     def _get(self, results_url):
#         """
#         Collect the results using GET. May have to try multiple times.
#
#         :return:
#         """
#         cmd = """curl {}""".format(results_url)
#
#         response = ast.literal_eval(self._run(cmd))
#         if response["status_code"] == 400:
#             print("status code: {}".format(response['status_code']))
#             raise RuntimeError()
#
#         status_code = int(response["status_code"])
#         while status_code > 200:
#             print("status code: {}".format(response['status_code']))
#             time.sleep(10)
#             response = ast.literal_eval(self._run(cmd))
#             status_code = int(response["status_code"])
#
#         print(type(response["status_code"]))
#         return response
#
#
# def main():
#     protoss = Protoss()
#     result = protoss.add_hydrogens(pdb_code="3cqw")
#
#     print(result.protein)
#
# if __name__ == "__main__":
#     main()