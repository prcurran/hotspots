"""
python wrapper for the rcsb PDB

Thank-you pypdb.py for a great starting point [ref]

- aiming to open up all search combinations
- composite queries
"""
from __future__ import print_function, division

import os
import sys
import shutil

import xmltodict
from hotspots.data import common_solvents
from tqdm import tqdm

if sys.version_info.major == 2:
    from urllib2 import Request as request
    from urllib2 import urlopen as urlopen

else:
    from urllib.request import Request as request
    from urllib.request import urlopen as urlopen


base_url = 'http://www.rcsb.org/pdb/rest'


class Helper(object):
    """
    generic functionality to aid RESTful python splender
    """

    @staticmethod
    def url_request(extension, data=None):
        """"""
        if data:
            url = request("{}/{}".format(base_url, extension), data=data)
        else:
            url = request("{}/{}".format(base_url, extension))

        return urlopen(url)


class Query(Helper):
    """
    constructs a query
    """

    class _QueryTerm(object):
        """

        """
        def __init__(self, query_type, query_parameters, conjunction):
            self.query_parameters = query_parameters
            self.conjunction = conjunction
            self.prefix = "org.pdb.query.simple."

            self.query_type = self.prefix + query_type
            self.query_parameters["queryType"] = self.query_type

            self._generate_description()
            self.term = {'orgPdbQuery': self.query_parameters}

        def _generate_description(self):
            """

            :return:
            """
            self.query_parameters["description"] = 'No description'

    def __init__(self):
        self.query_terms = []

    @property
    def query(self):
        """
        Assembles the QueryTerms added into correct dictionary format
        :return:
        """
        query_term_dict = {"queryRefinement": []}
        for i, query_term in enumerate(self.query_terms):

            if i > 0:
                try:
                    query_term.term.update({"conjunctionType": query_term.conjunction})
                except:
                    raise ValueError("For Queries with > 1 QueryTerm, conjuction can not be NoneType")

            query_term.term.update({"queryRefinementLevel": str(i)})
            query_term_dict["queryRefinement"].append(query_term.term)

        return {"orgPdbCompositeQuery": query_term_dict}

    def add_term(self, query_type, query_parameters, conjunction=None):
        """

        :param query_type:
        :param query_parameters:
        :return:
        """
        self.query_terms.append(Query._QueryTerm(query_type, query_parameters, conjunction))


class PDB(object):
    """
    this is a toolkit to access the PDB's RESTful web service.
    """
    @staticmethod
    def search(query=None):
        """
        search the pdb, if query is None the whole PDB will be returned
        :param query:
        :return:
        """
        print("Searching PDB ...")
        if query:
            query_text = xmltodict.unparse(query, pretty=False)
            query_text = query_text.encode()
            f = Helper.url_request(extension='/search', data=query_text)

            return [PDBResult(identifier=code.split(":")[0]) for code in tqdm(f.readlines())]

        else:
            f = Helper.url_request(extension='/getCurrent')
            f.read()


class _Ligand(Helper):
    """
    a class to handle PDB ligands (keep at string level)
    """

    def __init__(self, ligand_info):
        """

        :param ligand_str:
        """
        self.ligand_info = ligand_info

    @staticmethod
    def from_chemicalid(chemicalid):
        """

        :param chemicalid:
        :return:
        """

        extension = "describeHet?chemicalID={}".format(chemicalid)
        data = xmltodict.parse(Helper.url_request(extension=extension),
                               attr_prefix='',
                               dict_constructor=dict)

        return _Ligand(ligand_info=data.values()[0].values()[0]['ligand'])


    def search(self, search_type, tanimoto=0.7):
        """
        will do a chemical structure search using _Ligands smile string.
        Search Types:
            - 'exact'
            - 'substructure'
            - 'superstructure'
            - 'similarity'

        :param search_type: type of search to be executed
        :param tanimoto: float, tanimoto cooefficent for similarity search. 0 = most dissimilar, 1 = identifical
        :return:
        """
        supported = ["exact", "substructure", "superstructure", "similarity"]
        assert search_type in supported

        extension = "smilesQuery?smiles={0}&search_type={1}".format(self.smiles, search_type)
        if search_type is 'similarity':
            extension += "&similarity={}".format(tanimoto)

        data = xmltodict.parse(self.url_request(extension=extension),
                               attr_prefix='',
                               dict_constructor=dict)

        for key, value in data.values()[0].items():
            if key == "ligandInfo":
                return [_Ligand(ligand_info=v) for v in value.values()[0]]

    @property
    def structure_id(self):
        return self.ligand_info["structureId"]

    @property
    def smiles(self):
        return self.ligand_info["smiles"]

    @property
    def InChIKey(self):
        return self.ligand_info["InChIKey"]

    @property
    def molecular_weight(self):
        return self.ligand_info["molecularWeight"]

    @property
    def chemical_id(self):
        return self.ligand_info["chemicalID"]

    @property
    def chemical_name(self):
        return self.ligand_info["chemicalName"]

    @property
    def formula(self):
        return self.ligand_info["formula"]

    @property
    def chemical_type(self):
        return self.ligand_info["type"]


class _Protein(object):
    """
    """
    class _SubUnit(object):
        def __init__(self, data):
            self.data = data

        @property
        def entity_number(self):
            return self.data["entityNr"]

        @property
        def chain(self):
            return self.data["chain"]["id"]

        @property
        def taxonomy(self):
            t = self.data["Taxonomy"]
            if type(t) is list:
                return [a["name"] for a in t]
            else:
                return t["name"]

        @property
        def name(self):
            try:
                if type(self.data["macroMolecule"]) is list:
                    print(self.data["macroMolecule"])
                    upr = [a["name"] for a in self.data["macroMolecule"] if a["accession"]["id"] != "P00720"]
                    return upr[0]
                else:
                    return self.data["macroMolecule"]["name"]

            except KeyError:
                return None

        @property
        def accession_id(self):
            try:
                if type(self.data["macroMolecule"]) is list:
                    print(self.data["macroMolecule"])
                    upr = [a["accession"]["id"] for a in self.data["macroMolecule"] if a["accession"]["id"] != "P00720"]
                    return upr[0]
                else:
                    return self.data["macroMolecule"]["accession"]["id"]

            except KeyError:
                return None


    def __init__(self, data):
        """
        :param ligand_str:
        """
        self.data = data
        self.sub_unit_data = []
        for data_item in self.data.values():
            for i in data_item.values():
                if type(i) is dict:
                    self.sub_unit_data.append(i)
                elif type(i) is unicode:
                    self.structure_id = i
                elif type(i) is list:
                    self.sub_unit_data.append(i[0])
                else:
                    raise TypeError("{}".format(i))

    def sub_unit_count(self):
        return len(self.sub_unit_data)

    @property
    def sub_units(self):
        return [self._SubUnit(sud) for sud in (self.sub_unit_data)]


class _EntryDescription(object):
    """"""
    def __init__(self, data):
        """

        :param ligand_str:
        """
        self.data = data

    @property
    def status(self):
        return self.data["status"]

    @property
    def residue_count(self):
        return self.data["nr_residues"]

    @property
    def structure_authors(self):
        return self.data["structure_authors"]

    @property
    def citation_authors(self):
        return self.data["citation_authors"]

    @property
    def resolution(self):
        return self.data["resolution"]

    @property
    def title(self):
        return self.data["title"]

    @property
    def release_date(self):
        return self.data["release_date"]

    @property
    def experiment_method(self):
        return self.data["expMethod"]

    @property
    def pubmedId(self):
        return self.data["pubmedId"]

    @property
    def keywords(self):
        return self.data["keywords"]

    @property
    def atom_count(self):
        return self.data["nr_atoms"]

    @property
    def entities(self):
        return self.data["nr_entities"]

    @property
    def last_modification_date(self):
        return self.data["last_modification_date"]

    @property
    def deposition_date(self):
        return self.data["deposition_date"]

    @property
    def related_structures(self):
        return [d["pdbId"] for d in self.data["relatedPDB"]]


class PDBResult(object):
    """
    class to handle operations on a single PDB code
    """

    def __init__(self, identifier):
        """

        :param identifier:
        """
        self.identifier = identifier

        self._ligands = self.get_ligands()

    def _raw_properties(self, info_type='ligand'):
        """

        :return:
        """
        info_type_dict = {'describe_pdb': '/describePDB?structureId=',
                          'describe_mol': '/describeMol?structureId=',
                          'ligand': '/ligandInfo?structureId=',
                          'pfam': '/hmmer?structureId=',
                          'general': '/getEntityInfo?structureId='}

        url = request(base_url + info_type_dict[info_type] + self.identifier)
        return urlopen(url)

    def custom(self, field="classification"):
        """
        return custom information
        link: https://www.rcsb.org/pdb/results/reportField.do


        :param field:
        :return:
        """
        extension = "/customReport.xml?pdbids={}&customReportColumns={}&service=wsfile&format=xml".format(self.identifier, field)
        url = request(base_url + extension)

        x = (urlopen(url))
        data = xmltodict.parse(x.read(), attr_prefix='', dict_constructor=dict)
        return data["dataset"]["record"]["dimStructure.{}".format(field)]

    @property
    def descriptions(self):
        """returns description"""
        data = xmltodict.parse(self._raw_properties(info_type='describe_pdb').read(),
                               attr_prefix='',
                               dict_constructor=dict)
        for pdb_description in data.values():
            data_item = pdb_description.values()[0]

        return _EntryDescription(data_item)

    @property
    def protein(self):
        """"""
        data = xmltodict.parse(self._raw_properties(info_type='describe_mol').read(),
                               attr_prefix='',
                               dict_constructor=dict)
        return _Protein(data.values()[0])

    @property
    def ligands(self):
        """"""
        return self._ligands

    def get_ligands(self):
        """returns Ligand"""
        data = xmltodict.parse(self._raw_properties(info_type='ligand').read(),
                               attr_prefix='',
                               dict_constructor=dict)

        ligs = []
        for structureId in data.values():
            for ligandInfo in structureId.items():
                for ligand in ligandInfo:
                    if type(ligand) is dict:
                        for data_item in ligand.values():
                            if type(data_item) == list:
                                for d in data_item:
                                    ligs.append(_Ligand(d))
                            else:
                                ligs.append(_Ligand(data_item))
        return ligs

    @property
    def filtered_ligands(self):
        """"""
        cs = common_solvents()
        return [lig for lig in self._ligands if lig.chemical_id not in cs and
                                                float(lig.molecular_weight) > 120]

    def download(self, out_dir, compressed=False, biological_assembly=False):
        """
        only PDB for now, I will put the rest in later #agile
        :param compressed:
        :param biological_assembly
        :return:
        """
        if biological_assembly:
            extension = ".pdb1"

        else:
            extension = ".pdb"

        url = "https://files.rcsb.org/download/{}{}".format(self.identifier, extension)
        pdbfile = urlopen(url).read()

        if compressed:
            out = os.path.join(out_dir, self.identifier)
            if not os.path.exists(out):
                os.mkdir(out)

            self.fname = os.path.join(out, self.identifier + extension)

            with open(self.fname, "w") as writer:
                writer.write(pdbfile)

            shutil.make_archive(out, 'zip', out)
            shutil.rmtree(out)

        else:
            self.fname = os.path.join(out_dir, self.identifier + extension)
            with open(self.fname, "w") as writer:
                writer.write(pdbfile)


