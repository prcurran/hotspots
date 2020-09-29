"""
The purpose of this module is to wrap around the Arpeggio algorithm and transform the data into
a format useful for the Fragment Hotspot Pharmacophore analysis.

- Arpeggio

"""
import logging
import operator
import os
import re
import sys
import tempfile
import traceback
from collections import OrderedDict
from functools import reduce
from shutil import copyfile

import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from ccdc.molecule import Coordinates
from ccdc.protein import Protein

from hotspots.wrapper_pdb import PDBResult

PDB_LINE_TEMPLATE = '{record: <6}{serial: >5} {atom_name: ^4}{altloc: ^1}{resname: ^3} {chain_id: ^1}{resnum: >4}{icode: ^1}   {x: >8.3f}{y: >8.3f}{z: >8.3f}{occ: >6.2f}{tfac: >6.2f}          {element: >2}{charge: >2}'


def clean(pdb_path, dry=None):
    """
    TODO: reimplement with CCDC API

    taken directly from Harry Jubb. THANKS!!
    https://github.com/harryjubb/pdbtools

    :param pdb_path:
    :return:
    """
    pdb_noext, pdb_ext = os.path.splitext(pdb_path)
    pdb_ext = pdb_ext.replace('.', '')

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(os.path.split(os.path.splitext(pdb_path)[0])[1], pdb_path)

    # REMOVE MULTIPLE MODELS
    # BY TAKING THE FIRST MODEL
    model = structure[0]
    output_label = "clean"

    # RAISE AN ERROR FOR TOO MANY ATOMS
    if len(list(model.get_atoms())) > 99999:
        try:
            raise ValueError('More than 99999 atoms in the PDB model!')
        except:
            traceback.print_exc(file=sys.stdout)
            exit(9)

    # DETERMINE POLYPEPTIDES AND CHAIN BREAKS
    ppb = PPBuilder()
    polypeptides = ppb.build_peptides(model, aa_only=False)

    # MAKE DATA STRUCTURES FOR CHAIN POLYPEPTIDES
    chain_ids = set([x.id for x in model.child_list])
    chain_pieces = OrderedDict()
    chain_polypeptides = OrderedDict()
    chain_break_residues = OrderedDict()
    chain_sequences = OrderedDict()

    for chain_id in chain_ids:
        chain_pieces[chain_id] = 0
        chain_break_residues[chain_id] = []
        chain_polypeptides[chain_id] = []

    # GET ALL POLYPEPTIDE RESIDUES IN THE MODEL
    polypeptide_residues = []

    for pp in polypeptides:
        for res in pp:
            polypeptide_residues.append(res)

    # GET THE CHAIN_ID(S) ASSOCIATED WITH EACH POLYPEPTIDE
    polypeptide_chain_id_sets = [set([k.get_parent().id for k in x]) for x in polypeptides]

    for e, polypeptide_chain_id_set in enumerate(polypeptide_chain_id_sets):

        # WARN IF NOT JUST ONE CHAIN ID ASSOCIATED WITH THE POLYPEPTIDE
        if len(polypeptide_chain_id_set) != 1:
            logging.warn('A polypeptide had {} chains associated with it: {}'.format(len(polypeptide_chain_id_set),
                                                                                     polypeptide_chain_id_set))

        for polypeptide_chain_id in polypeptide_chain_id_set:
            chain_pieces[polypeptide_chain_id] = chain_pieces[polypeptide_chain_id] + 1

            # ADD FIRST AND LAST RESIDUE TO THE CHAIN BREAK RESIDUES (POLYPEPTIDE TERMINAL RESIDUES)
            chain_break_residues[polypeptide_chain_id] = chain_break_residues[polypeptide_chain_id] + [
                polypeptides[e][0], polypeptides[e][-1]]
            chain_polypeptides[polypeptide_chain_id] = chain_polypeptides[polypeptide_chain_id] + [polypeptides[e]]

    # POP OUT THE FIRST AND LAST RESIDUES FROM THE CHAIN BREAK RESIDUES
    # TO REMOVE THE GENUINE TERMINI
    for chain_id in chain_break_residues:
        chain_break_residues[chain_id] = chain_break_residues[chain_id][1:-1]

    all_chain_break_residues = reduce(operator.add, chain_break_residues.values())

    # MAKE THE CHAIN SEQUENCES FROM THE CHAIN POLYPEPTIDE PIECES
    for chain_id in chain_polypeptides:

        pp_seqs = [str(x.get_sequence()) for x in chain_polypeptides[chain_id]]

        if pp_seqs:
            chain_sequences[chain_id] = reduce(operator.add, pp_seqs)

    # WRITE OUT CLEANED PDB
    # MANY OF THE ISSUES ARE SOLVED DURING THE WRITING OUT
    print(f"{pdb_noext}_clean.{pdb_ext}")
    with open(f"{pdb_noext}_clean.{pdb_ext}", 'w') as fo:

        atom_serial = 1

        for residue in model.get_residues():
            record = 'ATOM'

            # REMOVE WATERS IF FLAG SET
            if dry:
                if residue.get_full_id()[3][0] == 'W':
                    continue

            # SET HETATM RECORD IF IT WAS ORIGINALLY A HETATM OR WATER
            if residue.get_full_id()[3][0] == 'W' or residue.get_full_id()[3][0].startswith('H_'):
                record = 'HETATM'

            # SET ATOM RECORD IF THE RESIDUE IS IN A POLYPEPETIDE
            if residue in polypeptide_residues:
                record = 'ATOM'

            # LOOP THROUGH ATOMS TO OUTPUT
            for atom in residue.child_list:

                # DEAL WITH DISORDERED ATOMS
                if atom.is_disordered():
                    atom = atom.disordered_get()

                # ALWAYS KEEP HYDROGEN

                # CONVERT SELENOMETHIONINES TO METHIONINES
                if residue in polypeptide_residues and (residue.resname == 'MSE' or residue.resname == 'MET'):

                    residue.resname = 'MET'

                    if atom.name == 'SE' and atom.element == 'SE':
                        atom.name = 'SD'
                        atom.element = 'S'

                # FIX ATOM NAME BUG
                if len(atom.name) == 3:
                    atom.name = ' ' + atom.name

                # PDB OUTPUT
                # ATOM SERIALS ARE RENUMBERED FROM 1
                # ALTLOCS ARE ALWAYS BLANK
                # CHARGES ARE ALWAYS BLANK(?)
                # OCCUPANCIES ARE ALWAYS 1.00
                output_line = PDB_LINE_TEMPLATE.format(record=record,
                                                       serial=atom_serial,
                                                       atom_name=atom.name,
                                                       altloc=' ',
                                                       resname=residue.resname,
                                                       chain_id=residue.get_parent().id,
                                                       resnum=residue.get_id()[1],
                                                       icode=residue.get_id()[2],
                                                       x=float(atom.coord[0]),
                                                       y=float(atom.coord[1]),
                                                       z=float(atom.coord[2]),
                                                       occ=1.00,
                                                       tfac=atom.bfactor,
                                                       element=atom.element,
                                                       charge='')

                fo.write('{}\n'.format(output_line))

                atom_serial += 1

                # RAISE AN ERROR IF WE'VE NOW GOT TOO MANY ATOMS
                if (atom_serial - 1) > 99999:

                    try:
                        raise ValueError('More than 99999 atoms in the PDB when renumbered!')
                    except:
                        traceback.print_exc(file=sys.stdout)
                        exit(9)


def _get_protein(pdb_code, tmpdir, protein_path=None):
    """
    Fetch the protein of interest

    NB: 13-May-2020: There seems to be some funny behaviour from the CCDC pdbwriter
    Therefore, AVOID writing using MoleculeWriter. (for now)

    :param str pdb_code: 4 character PDB code
    :param str tmpdir: path to temporary directory
    :param str protein_path: path to protein file
    :return: `ccdc.protein.Protein`
    """
    ext = "pdb"
    if protein_path:
        # copy protein to the results 'tmp' directory
        ext = protein_path.split(".")[1]
        copyfile(protein_path, os.path.join(tmpdir, f"{pdb_code}.{ext}"))

    else:
        # otherwise download from the PDB
        PDBResult(identifier=pdb_code).download(out_dir=tmpdir)

    # clean(os.path.join(tmpdir, f"{pdb_code}.pdb"), dry)  # clean suffix added
    return Protein.from_file(os.path.join(tmpdir, f'{pdb_code}.{ext}'))


def _get_ligand(prot, hetid, chain):
    """
    Fetch the ligand of interest

    :param `ccdc.protein.Protein` prot: protein containing the ligand.
    :param str hetid: 3 character ligand identifier
    :param str chain: chain identifier
    :return: `ccdc.molecule.Molecule`
    """
    prot.detect_ligand_bonds()
    return [l for l in prot.ligands
            if l.identifier.split(":")[0] == chain and l.identifier.split(":")[1][:3] == hetid][0]


class Arpeggio:
    """
    links

    """

    # TODO: Using arpeggio directly is sufficient for this work.
    #       However, reimplementation for smoother integration
    #       with ccdc.pharmacophore module would be nice.
    def __init__(self, pdb_code, hetid, chain="A", protein_path=None, tmpdir=None):
        if not tmpdir:
            self.tmpdir = tempfile.mkdtemp()
        else:
            self.tmpdir = tmpdir

        # 1) input data
        self.pdb_code = pdb_code
        self.hetid = hetid
        self.chain = chain

        self.protein_path = protein_path

        self.protein = _get_protein(pdb_code, self.tmpdir, protein_path)

        self.ligand = _get_ligand(self.protein, hetid, chain)

        self.ligand_label_dic = {a.label: i for i, a in enumerate(self.ligand.atoms)}
        self.prot_chain = [c for c in self.protein.chains if c.identifier == chain][0]

        ##################################################################################
        # the residue index and residue number are not equal
        # These dictionaries map the ccdc object index and the residue number
        self.protein_chain_index_dic = {r.identifier.split(":")[1][3:]: i
                                        for i, r in enumerate(self.prot_chain.residues)}
        self.protein_water_index_dic = {w.identifier.split(":")[1][3:]: j
                                        for j, w in enumerate(self.protein.waters)}
        self.protein_ligand_index_dic = {l.identifier.split(":")[1][3:]: k
                                        for k, l in enumerate(self.protein.ligands)}
        self.protein_metal_index_dic = {m.label: h
                                        for h, m in enumerate(self.protein.metals)}
        ##################################################################################

        self.ligand_index = int(self.ligand.identifier.split(":")[1][3:])

        self.arpeggio_dic = {'hbond acceptor': "acceptor",
                             'weak hbond acceptor': "weak_acceptor",
                             'pos ionisable': "positive",
                             'neg ionisable': "negative",
                             'aromatic': "apolar",
                             'carbonyl oxygen': None,
                             'carbonyl carbon': None,
                             'hydrophobe': "apolar",
                             'hbond donor': "donor",
                             'xbond acceptor': None,
                             'weak hbond donor': "weak_donor"}

    def run(self):
        """
        Run Arpeggio from the docker image
        :return: None
        """
        cmd = f"""
        docker run --rm -v "{self.tmpdir}":/run -u `id -u`:`id -g` harryjubb/arpeggio python arpeggio.py /run/{self.pdb_code}.pdb -s RESNAME:{self.hetid} -v
        """
        os.system(cmd)

        # After run there should be files

    def _get_summary(self):
        """
        Read and process the *.sift

        :return: `pandas.DataFrame`
        """
        headers = ["atom",
                   "clash",
                   "covalent",
                   "vdwclash",
                   "vdw",
                   "proximal",
                   "hbond",
                   "weakhbond",
                   "halogenbond",
                   "ionic",
                   "metalic",
                   "aromatic",
                   "hydrophobic",
                   "carbonyl",
                   "polar",
                   "weakpolar"
                   ]

        summary = pd.read_csv(os.path.join(self.tmpdir, f"{self.pdb_code}.sift"),
                              sep='\t',
                              names=headers,
                              index_col=False)
        summary["resid"] = [int(str(row.atom).split("/")[1]) for index, row in summary.iterrows()]
        summary["atomid"] = [str(row.atom).split("/")[2] for index, row in summary.iterrows()]
        summary = summary.loc[summary.resid == self.ligand_index]
        summary.reset_index()
        return summary

    def _get_contacts(self):
        """
        Read and process the *.contacts

        :return: `pandas.DataFrame`
        """
        headers = ["atom1",
                   "atom2",
                   "clash",
                   "covalent",
                   "vdwclash",
                   "vdw",
                   "proximal",
                   "hbond",
                   "weakhbond",
                   "halogenbond",
                   "ionic",
                   "metalic",
                   "aromatic",
                   "hydrophobic",
                   "carbonyl",
                   "polar",
                   "weakpolar",
                   "relationship"]

        contacts = pd.read_csv(os.path.join(self.tmpdir, f"{self.pdb_code}.contacts"),
                               sep='\t',
                               names=headers,
                               index_col=False)

        # str to convert pandas.Series to str, int to convert str to int for self.ligand_index evaluation.
        res_1 = [int(str(row["atom1"]).split("/")[1]) for index, row in contacts.iterrows()]
        res_2 = [int(str(row["atom2"]).split("/")[1]) for index, row in contacts.iterrows()]
        contacts["residue1"] = res_1
        contacts["residue2"] = res_2

        contacts = contacts.loc[(contacts.residue1 == self.ligand_index) | (contacts.residue2 == self.ligand_index)]
        contacts = contacts.reset_index()

        atomid = []
        for index, row in contacts.iterrows():
            if row.residue1 == self.ligand_index:
                lig_atom = str(row.atom1).split("/")[2]
            elif row.residue2 == self.ligand_index:
                lig_atom = str(row.atom2).split("/")[2]
            else:
                lig_atom = "NaN"
            atomid.append(lig_atom)

        contacts["atomid"] = atomid

        return contacts

    def _get_atom_types(self):
        """
        Read and process the *.atomtypes

        TODO: Atom-typing is done through SMARTS matches. The SMARTS definitions
              should be the same as the CCDC pharmacophore module
        :return: `pandas.DataFrame`
        """
        atom_types = pd.read_csv(os.path.join(self.tmpdir, f"{self.pdb_code}.atomtypes"),
                                 sep='\t',
                                 names=["atom", "atomtype"],
                                 index_col=False)

        atom_types.atomtype = [re.findall(r"'(.*?)'", row.atomtype) for index, row in atom_types.iterrows()]
        return atom_types

    def _get_rings(self):
        """
        Read and process the *.rings

        :return: `pandas.DataFrame`
        """
        rings = pd.read_csv(os.path.join(self.tmpdir, f"{self.pdb_code}.rings"),
                            sep='\t',
                            names=["ringid", "resid", "centroid"],
                            index_col=False)
        rings.centroid = [[float(a) for a in str(row.centroid).strip("][").split(",")] for index, row in
                          rings.iterrows()]

        rings["resnum"] = [int(str(row["resid"]).split("/")[1]) for index, row in rings.iterrows()]
        return rings

    def _get_ring_to_ring(self):
        """
        Read and process the *.ri

        :return: `pandas.DataFrame`
        """

        ring_to_ring = pd.read_csv(os.path.join(self.tmpdir, f"{self.pdb_code}.ri"),
                                   sep='\t',
                                   names=["ringid1", "resid1", "centroid1",
                                          "ringid2", "resid2", "centroid2",
                                          "interactiontypes", "relationshipres", "relationship"],
                                   index_col=False)
        ring_to_ring.centroid1 = [[float(a) for a in str(row.centroid1).strip("][").split(",")]
                                  for index, row in ring_to_ring.iterrows()]
        ring_to_ring.centroid2 = [[float(a) for a in str(row.centroid2).strip("][").split(",")]
                                  for index, row in ring_to_ring.iterrows()]
        return ring_to_ring

    def _get_atom_to_ring(self):
        """
        Read and process the *.ari

        :return: `pandas.DataFrame`
        """
        atom_to_ring = pd.read_csv(os.path.join(self.tmpdir, f"{self.pdb_code}.ari"),
                                   sep='\t',
                                   names=["atom", "ringid", "resid", "centroid", "interactiontypes", "atype",
                                          "relationship"],
                                   index_col=False)
        atom_to_ring["atomid"] = [str(row.atom).split("/")[2] for index, row in atom_to_ring.iterrows()]
        atom_to_ring["resnum"] = [int(str(row.atom).split("/")[1]) for index, row in atom_to_ring.iterrows()]
        atom_to_ring.centroid = [str(row.centroid).strip("][").split(",") for index, row in atom_to_ring.iterrows()]
        atom_to_ring.interactiontypes = [re.findall(r"'(.*?)'", row.interactiontypes) for index, row in
                                         atom_to_ring.iterrows()]
        return atom_to_ring

    def _get_partner_atom(self, atm1, atm2, atm):
        """
        Given a contact pair, "protein" atom interacting with the ligand is returned

        "protein" includes waters, other ligands and metals.
        :param atm1:
        :param atm2:
        :param atm:
        :return:
        """
        # we don't know the order of atoms in output, so remove the ligand atom by id

        prot_hbond_atm = [atm1, atm2]
        prot_hbond_atm.remove(f"{self.chain}/{self.ligand_index}/{atm.label}")
        # split
        c, p_residue_index, p_atm_label = prot_hbond_atm[0].split("/")

        # identify the protein atom type
        p_residue_index = str(p_residue_index)
        if p_residue_index in self.protein_chain_index_dic:
            print(self.protein_chain_index_dic[p_residue_index])
            p_residue_obj = self.prot_chain.residues[self.protein_chain_index_dic[p_residue_index]]
            p_atom_obj = [a for a in p_residue_obj.atoms if a.label == p_atm_label]

        elif p_residue_index in self.protein_water_index_dic:
            p_residue_obj = self.protein.waters[self.protein_water_index_dic[p_residue_index]]
            p_atom_obj = [a for a in p_residue_obj.atoms if a.atomic_symbol is not "H"]

        elif p_residue_index in self.protein_ligand_index_dic:
            p_residue_obj = self.protein.ligands[self.protein_ligand_index_dic[p_residue_index]]
            p_atom_obj = [a for a in p_residue_obj.atoms if a.label == p_atm_label]

        elif p_residue_index in self.protein_metal_index_dic:
            p_residue_obj = self.protein.metals[self.protein_metal_index_dic[p_residue_index]]
            p_atom_obj = [a for a in p_residue_obj.atoms if a.label == p_atm_label]

        else:
            p_atom_obj = None

        return prot_hbond_atm, p_atom_obj

    def get_hbonds(self, atm, relationship=None):
        """
        Get Hbond pairing, atom positions and objects

        :param `ccdc.molecule.Atom` atm: a CCDC Atom
        :return: list of `hotspots.wrapper_arpeggio.Interaction`
        """
        features = []
        contacts = self._get_contacts()
        atom_types = self._get_atom_types()

        # atoms of interest
        if relationship:
            atm_hbonds = contacts.loc[(contacts.atomid == atm.label) &
                                      (contacts.hbond == 1) &
                                      (contacts.relationship == relationship)]
        else:
            atm_hbonds = contacts.loc[(contacts.atomid == atm.label) & (contacts.hbond == 1)]

        if len(atm_hbonds) == 0:
            return features

        # grab the point `ccdc.molecule.Atom`
        lig_atm_type = atom_types.loc[atom_types["atom"] == f"{self.chain}/{self.ligand_index}/{atm.label}"]

        reclassified = [self.arpeggio_dic[i] for i in set(lig_atm_type.atomtype.values[0])
                        if i == "hbond donor"
                        or i == "hbond acceptor"]

        for index, row in atm_hbonds.iterrows():
            # grab the projected `ccdc.molecule.Atom`
            prot_hbond_atm, p_atom_obj = self._get_partner_atom(str(row.atom1), str(row.atom2), atm)

            # determine whether ligand atom is hbond is donor or acceptor
            if len(reclassified) == 1:
                con_type = reclassified[0]

            elif len(reclassified) == 2:
                # if the ligand atom is a doneptor check the protein atom
                if p_atom_obj[0].is_donor and not p_atom_obj[0].is_acceptor:
                    con_type = "acceptor"
                elif p_atom_obj[0].is_acceptor and not p_atom_obj[0].is_donor:
                    con_type = "donor"
                # elif angle is good == donor
                else:
                    con_type = "doneptor"
                    print("doneptor")
                    # TODO: at this point we would need to look at H atom placement
                    #       come back to this if needed
            else:
                raise TypeError

            feat = Interaction(point=atm.coordinates,
                               point_atom=atm,
                               point_identifier=f"{self.hetid}/{atm.label}",
                               point_type=con_type,
                               projected=p_atom_obj[0].coordinates,
                               projected_atom=p_atom_obj[0],
                               projected_identifier=prot_hbond_atm[0],
                               interaction_type="hbond")

            features.append(feat)
        return features

    def get_weak_hbonds(self, atm, relationship=None):
        """
        Get weakhbond pairing, atom positions and objects. Only returns weakhbond if atom is not participating in
        a 'proper' hbond.

        :param `ccdc.molecule.Atom` atm: a CCDC Atom
        :return: list of `hotspots.wrapper_arpeggio.Interaction`
        """
        features = []
        contacts = self._get_contacts()
        atom_types = self._get_atom_types()

        # NB: only interested in weak hbond contacts if a proper hbond is not present
        # obviously, this is not True for every use.
        if relationship:
            atm_weak_hbonds = contacts.loc[(contacts.atomid == atm.label) &
                                           (contacts.relationship == relationship) &
                                           ((contacts.hbond == 0) & (contacts.weakhbond == 1))]
        else:
            atm_weak_hbonds = contacts.loc[(contacts.atomid is atm.label) &
                                           ((contacts.hbond == 0) & (contacts.weakhbond == 1))]

        if len(atm_weak_hbonds) == 0:
            return features

        lig_atm_type = atom_types.loc[atom_types["atom"] == f"{self.chain}/{self.ligand_index}/{atm.label}"]

        # atom could still be type weak or strong donor / acceptor
        # FHM doesn't evaluate placement of weak acceptors / weak donors (perhaps remove them???)
        #           - weak donor / acceptor on the protein is fine
        #           - weak donor / acceptor on the ligand is iffy

        converted_lig_atm_type = [self.arpeggio_dic[i] for i in set(lig_atm_type.atomtype.values[0])
                                  if
                                  i == "weak hbond donor" or
                                  i == "weak hbond acceptor" or
                                  i == "hbond donor" or
                                  i == "hbond acceptor"]

        for index, row in atm_weak_hbonds.iterrows():
            prot_hbond_atm, p_atom_obj = self._get_partner_atom(str(row.atom1), str(row.atom2), atm)
            print(prot_hbond_atm, p_atom_obj)
            feat = Interaction(point=atm.coordinates,
                               point_atom=atm,
                               point_identifier=f"{self.hetid}/{atm.label}",
                               point_type=f"{converted_lig_atm_type[0]}",
                               projected=p_atom_obj[0].coordinates,
                               projected_atom=p_atom_obj[0],
                               projected_identifier=prot_hbond_atm[0],
                               interaction_type="weak_hbond")

            features.append(feat)
        return features

    def get_ligand_atom_to_ring_bonds(self, atm):
        """
        Get atom to ring pairing, atom positions and objects

        NB: contact types have been simplified, for a comparison to fragment hotspot maps which can
        not differentiate between aromatic contract types

        :param `ccdc.molecule.Atom` atm: a CCDC Atom
        :return: list of `hotspots.wrapper_arpeggio.Interaction`
        """
        features = []
        atom_to_ring = self._get_atom_to_ring()
        if len(atom_to_ring) == 0:
            return features

        atom_to_ring.to_csv(os.path.join(self.tmpdir, "test.csv"))
        atm_to_ring_contact = atom_to_ring.loc[(atom_to_ring.resnum == self.ligand_index) &
                                               (atom_to_ring.atomid == atm.label)]

        for index, row in atm_to_ring_contact.iterrows():
            # TODO: it would be nice to return the projected atoms of a ring not essential for this work though
            feat = Interaction(point=atm.coordinates,
                               point_atom=atm,
                               point_identifier=f"{self.hetid}/{atm.label}",
                               point_type="atom",
                               projected=Coordinates(x=row.centroid[0], y=row.centroid[1], z=row.centroid[2]),
                               projected_atom=None,
                               projected_identifier=f"{row.resid}ringid_{row.ringid}",
                               interaction_type="aromatic")

            try:
                res = self.protein_chain_index_dic[str(row.resnum)]
                feat.projected_residue = self.prot_chain.residues[res]
            except:
                # do nothing
                print(f"Residue: {row.resnum} is not in Chain {self.chain}")

            feat.feature_type_specific = str(row.interactiontypes)

            features.append(feat)
        return features

    def get_ligand_ring_to_atom_bonds(self, ringid, relationship=None):
        """
        Get ring to atom pairing, atom positions and objects

        NB: contact types have been simplified
        :param str ringid:
        :return: list of `hotspots.wrapper_arpeggio.Interaction`
        """
        features = []
        rings = self._get_rings()
        atom_to_ring = self._get_atom_to_ring()
        if len(atom_to_ring) == 0:
            return features

        centroid = rings.loc[rings.ringid == ringid].centroid.values[0]

        if relationship:
            ring_to_atom_contact = atom_to_ring.loc[(atom_to_ring.ringid == ringid) &
                                                    (atom_to_ring.relationship == relationship)]
        else:
            ring_to_atom_contact = atom_to_ring.loc[(atom_to_ring.ringid == ringid)]

        for index, row in ring_to_atom_contact.iterrows():
            chain, resnum, label = str(row.atom).split("/")
            try:
                res = self.prot_chain.residues[self.protein_chain_index_dic[resnum]]
                atm = [a for a in res.atoms if a.label == label][0]

            except:
                res = self.protein.waters[self.protein_water_index_dic[resnum]]
                atm = [a for a in res.atoms if a.atomic_symbol != "H"][0]

            feat = Interaction(point=Coordinates(x=centroid[0], y=centroid[1], z=centroid[2]),
                               point_atom=None,
                               point_identifier=f"{self.hetid}/ringid_{ringid}",
                               point_type="ring",
                               projected=atm.coordinates,
                               projected_atom=atm,
                               projected_identifier=str(row.atom),
                               interaction_type="aromatic")

            feat.feature_type_specific = str(row.interactiontypes)

            features.append(feat)
            print(features)
        return features

    def get_ring_to_ring_bonds(self, ringid):
        """
        Get ring to ring pairing, atom positions and objects

        NB: contact types have been simplified
        :param str ringid:
        :return: list of `hotspots.wrapper_arpeggio.Interaction`
        """
        features = []
        ring_to_ring = self._get_ring_to_ring()

        if len(ring_to_ring) == 0:
            return features

        ring_to_ring_contact = ring_to_ring.loc[(ring_to_ring.ringid1 == ringid) |
                                                (ring_to_ring.ringid2 == ringid)]

        pairings = []
        for index, row in ring_to_ring_contact.iterrows():
            if row.relationship == "INTER":
                # deduplicate, we only need the interaction one way around.
                # (Not sure if this is intentional behaviour?)
                new_pairing = {row.ringid1, row.ringid2}
                if any([new_pairing == a for a in pairings]):
                    continue
                else:
                    pairings.append({row.ringid1, row.ringid2})

                if row.ringid1 == ringid:
                    point, pointid = [row.centroid1, row.resid1]
                    projected, projectedid = [row.centroid2, row.resid2]
                    projectedring = row.ringid2
                else:                                                   # int(row.ringid2) == ringid
                    point, pointid = [row.centroid2, row.resid2]
                    projected, projectedid = [row.centroid1, row.resid1]
                    projectedring = row.ringid1

                feat = Interaction(point=Coordinates(x=point[0], y=point[1],  z=point[2]),
                                   point_atom=None,
                                   point_identifier=f"{self.hetid}/ringid_{ringid}",
                                   point_type="ring_to_ring",
                                   projected=Coordinates(x=projected[0], y=projected[1], z=projected[2]),
                                   projected_atom=None,
                                   projected_identifier=f"{projectedid}ringid_{projectedring}",
                                   interaction_type="aromatic")

                feat.feature_type_specific = str(row.interactiontypes)

                features.append(feat)
        return features

    def create_feature_list(self, relationship="INTER"):
        """

        :param str relationship: the relationship between the contact atoms. i.e INTER, INTRA, SELECTION_WATER ..
        :return:
        """
        atom_features = []

        # ligand atom centred contacts
        for atm in self.ligand.heavy_atoms:
            # for the purpose of comparison to fragment hotspot maps, we only let each atom take one
            # intermolecular interaction, this is obviously a simplification -- not yet sure if this is
            # an over-simplification
            feats = self.get_hbonds(atm, relationship=relationship)
            if len(feats) == 0:
                feats = self.get_weak_hbonds(atm, relationship=relationship)
                if len(feats) == 0:
                    feats = self.get_ligand_atom_to_ring_bonds(atm)
            atom_features.extend(feats)

        # ligand ring centred contacts
        ring_features = []
        rings = self._get_rings()
        ligand_rings = rings.loc[rings.resnum == self.ligand_index].ringid.values
        for ligand_ring in ligand_rings:
            # We will take all aromatic contact pairs, since on fragment hotspot maps the "apolar" propensity
            # is pretty diffused
            ring_features.extend(self.get_ligand_ring_to_atom_bonds(ligand_ring, relationship=relationship))
            ring_features.extend(self.get_ring_to_ring_bonds(ligand_ring))

        # deduplicate
        return atom_features, ring_features


class Interaction:
    """
    storage for the information that I am interested in.
    """

    def __init__(self, point, point_atom, point_identifier, point_type,
                 projected, projected_atom, projected_identifier,
                 interaction_type):
        self.point = point
        self.point_atom = point_atom
        self.point_identifier = point_identifier
        self.point_type = point_type

        self.projected = projected
        self.projected_atom = projected_atom
        self.projected_identifier = projected_identifier

        self.interation_type = interaction_type

    def __repr__(self):
        return f"{self.interation_type}: {self.point_identifier}-{self.projected_identifier}"
