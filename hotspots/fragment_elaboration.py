#ref https://iwatobipen.wordpress.com/2016/02/07/several-reaction-steps-handling-in-rdkit/

import os
from os import listdir

import numpy as np
from ccdc import io
from ccdc import molecule
from ccdc.docking import Docker
from ccdc.io import _CSDDatabaseLocator
from ccdc.protein import Protein
from ccdc.utilities import _private_importer
from ccdc_internal import molecule
from ccdc_internal.interaction import Grid
from rdkit import Chem
from rdkit.Chem import AllChem

with _private_importer():
    import ChemicalAnalysisLib
    import ConformerGeneratorLib

from fhm.fragment_hotspots import Hotspots

from collections import OrderedDict
from operator import itemgetter
import tempfile


class FragmentElaboration(object):
    def __init__(self, in_dir, charged= False, library = True):
        self.in_dir = in_dir

        self.fragment = os.path.join(self.in_dir, "fragment.mol2")
        if library == True:
            self.virtual_library = [m for m in io.MoleculeReader("virtual_library.mol2")]
        else:
            self.virtual_library = self.generate_library()


        self.protein = os.path.join(self.in_dir, "protein.pdb")
        self.charged = charged
        self.hotspot_result = self.generate_hotspot()
        self.bcv_result = self.generate_BCV()
        self.layer_dict = self.generate_layer()

        self.constraints = self.growing_constraints(self.layer_dict["2"])

        #self.add_ligands = os.path.join(self.in_dir, "decorated_fragments.mol2")
        #self.hit_list = self.run_docking()

        # {ccdc.Molecule: [percentage_overlap, hotspot_score]
        self.score_dict = {}

    def generate_library(self):
        print "Decorate fragment..."
        m = Chem.MolFromMol2File(self.fragment)

        #inputs (place saver)
        mol = Chem.MolToSmiles(m)
        terminal = Chem.SDMolSupplier("r1-20.sdf")
        spacer = Chem.SDMolSupplier("r2-20.sdf")

        vl = self.supplier(mol, terminal, spacer)
        mols = [self.from_smiles(v) for v in vl]

        with io.MoleculeWriter("decorated_fragments.mol2") as w:
            for m in mols:
                w.write(m)

        return mols

    def supplier(self, mol, terminal, spacer):

        term = ["({})".format(Chem.MolToSmiles(t)[3:]) for t in terminal]
        rxn = [self.aromaticsubstitution(t) for t in term]
        vl = [self.runreaction([mol], r) for r in rxn]
        virtual_library = [m for s in vl for m in s]
        return virtual_library

    def from_smiles(self, smiles, identifier=None, generate_initial_sites=True):
        if identifier is None:
            identifier = smiles

        if generate_initial_sites:
            parameter_files = _CSDDatabaseLocator.get_conformer_parameter_file_location()
            molmaker = ConformerGeneratorLib.MoleculeTo3D(parameter_files)
            mol = molecule.Molecule(identifier, molmaker.create_conformation(smiles))
        else:
            molmaker = ChemicalAnalysisLib.SMILESMoleculeMaker()
            mol = molecule.Molecule(identifier, _molecule=molmaker.siteless_atoms(smiles))
        return mol


    def aromaticsubstitution(self, smi):
        #example reaction
        smarts = "[cH&$(c(c)c):2]>>[c:2]{}".format(smi)
        rxn = AllChem.ReactionFromSmarts(smarts)
        return rxn

    def runreaction(self, list_smiles, rxn):
        mols = [Chem.MolFromSmiles(smi) for smi in list_smiles]
        products = set()
        for mol in mols:
            ps = rxn.RunReactants((mol,))
            for x in ps:
                products.add(Chem.MolToSmiles(x[0], isomericSmiles=True))
        return products



    def generate_hotspot(self):
        # generate HS from grd or from protein
        print "Run Fragment Hotspots..."
        h = Hotspots()

        inputs = [f for f in listdir(self.in_dir) if f.endswith(".grd") or f.endswith(".pdb")]

        if len([p for p in inputs if "protein.pdb" in p]) == 0:
            raise ValueError("No protein file found")

        if inputs > 1:
            if self.charged == False:
                return h.from_grid_dic(super_grids= {"apolar": Grid.from_file(os.path.join(self.in_dir, "apolar.grd")),
                                                     "donor": Grid.from_file(os.path.join(self.in_dir, "donor.grd")),
                                                     "acceptor": Grid.from_file(os.path.join(self.in_dir, "acceptor.grd"))
                                                     },
                                       prot= Protein.from_file(self.protein),

                                       )

            else:
                return h.from_grid_dic(super_grids={"apolar": Grid.from_file(os.path.join(self.in_dir, "apolar.grd")),
                                                    "donor": Grid.from_file(os.path.join(self.in_dir, "donor.grd")),
                                                    "acceptor": Grid.from_file(os.path.join(self.in_dir, "acceptor.grd")),
                                                    "negative": Grid.from_file(os.path.join(self.in_dir, "negative.grd")),
                                                    "positive": Grid.from_file(os.path.join(self.in_dir, "positive.grd"))
                                                    },
                                       prot=Protein.from_file(self.protein)
                                  )
        else:
            #superstar isn't currently working on windows
            return h.from_protein(prot=os.path.join(self.in_dir, "protein.pdb"), charged_probes=self.charged,
                                  experimental_buriedness=True)


    def generate_BCV(self):
        bcv = self.hotspot_result.best_continuous_volume(volume = 500)
        bcv._minimal_grid()

        if not os.path.exists(os.path.join(self.in_dir, "bcv")):
            os.mkdir(os.path.join(self.in_dir, "bcv"))

        os.chdir(os.path.join(self.in_dir, "bcv"))
        bcv.output_pymol_file()
        os.chdir(os.path.join(self.in_dir))

        return bcv

    def set_uniform_values(self, grd, threshold = 0, set = 1, mode= "less_than"):
        nx, ny, nz = grd.nsteps
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    if mode == "less_than":
                        if grd.value(i, j, k) < threshold:
                            grd.set_value(i, j, k, set)
                    else:
                        if grd.value(i, j, k) != threshold:
                            grd.set_value(i,j,k,set)

        return grd


    def _get_sphere_grid(self, template, molecule, size_scaling = 0.75):
        layer = template
        layer *= 0
        for a in molecule.atoms:
            layer.set_sphere(a.coordinates, a.vdw_radius * size_scaling, 1, scaling="none")
            layer = self.set_uniform_values(grd=layer, mode = 'other')
        return layer


    def _diff_to_map(self, diff, map, prot,):

        grd_dict = {}
        for probe, g in self.bcv_result.super_grids.items():
            mg = map[probe]
            overlap = (diff > 0) & (mg > 0)
            new = overlap * g
            grd_dict.update({probe: new})

        return Hotspots.Hotspot_results(grd_dict, prot, self.in_dir, ligsite_grid=None)


    def generate_layer(self, iterations = 2, thickness = 10):
        layer_dict = {}
        prot = Protein.from_file(self.protein)
        masked_grids = self.bcv_result._single_grid()

        for i in range(1,iterations +1):
            if i == 1:
                #initial layer
                difference_layer = self._get_sphere_grid(template=self.bcv_result.super_grids["apolar"].copy(),
                                                         molecule=io.MoleculeReader(self.fragment)[0])
                self.inner = difference_layer
                hr = self._diff_to_map(diff=difference_layer, map=masked_grids, prot=prot)

            else:
                self.outer = self.inner.dilate()
                for j in range(1,thickness):
                    self.outer = self.outer.dilate()

                difference_layer = self.outer - self.inner
                self.inner = self.outer

                hr = self._diff_to_map(diff=difference_layer, map=masked_grids, prot=prot)

            layer_dict.update({"{}".format(i): hr})

        return layer_dict


    def growing_constraints(self, hr):
        #return {polar: [constraint1, contstain2], apolar: [fitting_points]}

        constraints = {"polar": [], "apolar": ""}

        #apolar
        for probe, g in hr.super_grids.items():
            cutoff = np.percentile(a=self.get_scores(hr.super_grids[probe]), q=95)
            g = self.set_uniform_values(g, threshold = cutoff, set = 0)
            if probe == "apolar":
                hr.output_apolar_fitting_points(fname=os.path.join(self.in_dir, "hs_apolar_fitting_pts.mol2"),
                                                cutoff=cutoff)
                constraints["apolar"] = os.path.join(self.in_dir, "hs_apolar_fitting_pts.mol2")

            else:
                #prepare layer grid
                include = [isla for isla in g.islands(threshold=cutoff) if isla.count_grid() > 10]
                if len(include) == 0:
                    hr.super_grids[probe] = g * 0
                else:
                    hr.super_grids[probe] = Grid.super_grid(0, *include)

        os.chdir("Z:/denovo_tool/second_layer_edit")
        hr.output_pymol_file()
        os.chdir("Z:/denovo_tool")

        scored_atoms_dic = hr._score_protein_atoms()
        print scored_atoms_dic
        for score in sorted(scored_atoms_dic.keys(), reverse=True):
            atoms = []
            for atom, residue, atom_type in scored_atoms_dic[score]:
                atoms.append(atom)

            filtered_atoms = [a for a in atoms if a.atomic_number == 1 or a.is_acceptor]
            for a in filtered_atoms:
                print a.index, a.label, a.coordinates

            if len(filtered_atoms) > 0 and score > 10:
                print len(filtered_atoms), score



        #polar




        return "balH"


    def run_docking(self):
        #take virtual library and run docking
        print "Run GOLD docking ..."

        docker = Docker()
        settings = docker.settings

        self.start_ligand = io.MoleculeReader(os.path.join(self.in_dir, "fragment.mol2"))[0]
        tempd = tempfile.mkdtemp()

        settings.add_protein_file(os.path.abspath(self.protein))
        settings.binding_site = settings.BindingSiteFromPoint(settings.proteins[0], self.start_ligand.centre_of_geometry(), 10.0)
        settings.fitness_function = 'plp'
        settings.autoscale = 10.
        settings.output_directory = tempd
        #settings.output_directory = self.in_dir
        settings.output_file = "docked_ligands.mol2"
        settings.add_ligand_file(self.add_ligands, ndocks = 10)

        #setup constraints
        settings.add_constraint(settings.TemplateSimilarityConstraint(type= "all", template = self.start_ligand, weight= 150))

        settings.ProteinFileInfo().fitting_points_file("fname.mol2")
        #feed in layer2
        #self.hotspot_result.predict_protein_hbond_constraints(settings, weight = 100)

        results = docker.dock()

        #fragment = results.ligands[0]
        ligand_reader = results.ligands
        output_file = os.path.join(settings.output_directory, settings.output_file)
        docked_molecules = [m for m in io.MoleculeReader(os.path.join(tempd, output_file))]
        print docked_molecules

        return docked_molecules


    def get_scores(self, g):
        nx, ny, nz = g.nsteps
        return [g.value(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz) if g.value(i, j, k) != 0]


    def fudge_shift(self, a, b):
        return np.sqrt(((a.coordinates.x - b.coordinates.x)**2 +
                         (a.coordinates.y - b.coordinates.y)**2 +
                         (a.coordinates.z - b.coordinates.z)**2
                         ))



    def score_hitlist(self):
        #volume overlap between initial fragment and subtitution
        #self.start_ligand
        overlap_cutoff = 0.85
        lig = io.MoleculeReader(os.path.join(self.in_dir, "fragment.mol2"))[0]

        ref = self.bcv_result.super_grids["apolar"].copy()
        ref *= 0

        for a in lig.atoms:
            ref.set_sphere(a.coordinates, a.vdw_radius, 1)

        #placeholder atom rmsd needed
        ref_a = [at for at in lig.atoms if at.is_cyclic and at.is_donor and at.atomic_weight == 14.0067][0]
        for i, h in enumerate(self.hit_list):
            print i
            clean = self.bcv_result.super_grids["apolar"].copy()
            clean *= 0
            for b in h.atoms:
                clean.set_sphere(b.coordinates, b.vdw_radius, 1)

            overlap = (ref > 0) & (clean > 0)

            percentage_overlap = float(len(self.get_scores(overlap)))/float(len(self.get_scores(ref)))

            hit_b = [atm for atm in h.atoms if atm.is_cyclic and atm.is_donor and atm.atomic_weight == 14.0067][0]
            fudge_shift = self.fudge_shift(ref_a, hit_b)

            #print percentage_overlap, len(self.get_scores(ref)), len(self.get_scores(clean)), len(self.get_scores(overlap))

            if percentage_overlap > overlap_cutoff and fudge_shift < 2:

                hotspot_score = self.hotspot_result.score_ligand(h)
                self.score_dict.update({h: hotspot_score})


        d = OrderedDict(sorted(self.score_dict.items(), key= itemgetter(1), reverse=True))

        return d

def main():
    f = FragmentElaboration(in_dir = "Z:\denovo_tool")

    # if not os.path.exists(os.path.join(f.in_dir, "first_layer")):
    #     os.mkdir(os.path.join(f.in_dir, "first_layer"))
    #
    # os.chdir(os.path.join(f.in_dir, "first_layer"))
    # f.layer_dict["1"].output_pymol_file()
    # os.chdir(os.path.join(f.in_dir))
    #
    # if not os.path.exists(os.path.join(f.in_dir, "second_layer")):
    #     os.mkdir(os.path.join(f.in_dir, "second_layer"))
    #
    # os.chdir(os.path.join(f.in_dir, "second_layer"))
    # f.layer_dict["2"].output_pymol_file()
    # os.chdir(os.path.join(f.in_dir))

    # with io.MoleculeWriter("not_scored_hitlist.mol2") as w:
    #     for k in f.hit_list:
    #         w.write(k)
    #
    #
    # d = f.score_hitlist()
    # print d
    #
    # with io.MoleculeWriter("scored_hitlist.mol2") as w:
    #     for i, m in enumerate(d.keys()):
    #         print i, d[m]
    #         w.write(m)





if __name__ == "__main__":
    main()
