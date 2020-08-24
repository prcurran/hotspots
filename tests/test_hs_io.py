import unittest
from ccdc.protein import Protein
from ccdc.io import MoleculeReader
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.result import Results
from hotspots.grid_extension import Grid
from hotspots.calculation import Runner


from hotspots.wrapper_pymol import PyMOLFile, PyMOLCommands

from pprint import pprint
import numpy as np
import os
import tempfile
import zipfile


class TestHotspotReader(unittest.TestCase):
    def test_read(self):
        path = "testdata/hs_io/minimal_multi_all_grids/out.zip"
        with HotspotReader(path=path) as r:
            hr = r.read(identifier="hotspot-1")

        self.assertIsInstance(hr, Results)

        with HotspotReader(path=path) as r:
            hr = r.read()

        self.assertIsInstance(hr, list)


class TestRun():
    def test_generate_real(self):
        runner = Runner()
        hr = runner.from_pdb(pdb_code="2vta", buriedness_method='ghecom')
        settings = HotspotWriter.Settings()
        settings.output_superstar = True
        parent = "testdata/2vta"
        with HotspotWriter(parent) as w:
            w.write(hr)

    def generate_fake(self, buriedness=False, weighted=False, superstar=True):
        """
        create a small set of grids for testing

        :param buriedness:
        :param weighted:
        :param superstar:
        :return:
        """

        def populate_grid(template, num_spheres, radius=1, value=8, scaling='linear'):
            h = template.copy_and_clear()
            for i in range(1, num_spheres):
                x, y, z = [np.random.randint(low=2, high=ax - 2, size=1) for ax in h.nsteps]

                h.set_sphere(point=h.indices_to_point(x, y, z),
                             radius=radius,
                             value=value,
                             scaling=scaling)

            return h

        protein = Protein.from_file("testdata/6y2g_A/binding_site.pdb")
        mol = MoleculeReader("testdata/6y2g_A/A_mol.mol2")[0]
        g = Grid.initalise_grid([a.coordinates for a in mol.atoms])

        if buriedness:
            buriedness_grid = Grid.from_molecule(mol)
        else:
            buriedness_grid = None

        interactions = ["apolar", "donor", "acceptor"]

        super_grids = {p: populate_grid(template=g, num_spheres=3) for p in interactions}

        if superstar:
            superstar_grids = {p: populate_grid(template=g, num_spheres=3) for p in interactions}
        else:
            superstar_grids = None

        if weighted:
            weighted_superstar_grids = {p: populate_grid(template=g, num_spheres=3) for p in interactions}
        else:
            weighted_superstar_grids = None

        return Results(super_grids=super_grids,
                       protein=protein,
                       buriedness=buriedness_grid,
                       superstar=superstar_grids,
                       weighted_superstar=weighted_superstar_grids)

    def test_write_real_single(self):
        base = "testdata/1hcl"
        interactions = ["donor", "acceptor", "apolar"]
        super_grids = {p: Grid.from_file(os.path.join(base, f"{p}.grd")) for p in interactions}
        superstar_grids = {p: Grid.from_file(os.path.join(base, f"superstar_{p}.grd")) for p in interactions}
        buriedness = Grid.from_file(os.path.join(base, "buriedness.grd"))
        prot = Protein.from_file(os.path.join(base, "protein.pdb"))

        hr = Results(super_grids=super_grids,
                     protein=prot,
                     buriedness=buriedness,
                     superstar=superstar_grids)

        settings = HotspotWriter.Settings()
        settings.output_superstar = True
        with HotspotWriter("testdata/hs_io/minimal_all_grids_real", settings=settings) as w:
            w.write(hr)

    def test_write_fake_single(self):
        a = self.generate_fake(buriedness=True, superstar=True)
        settings = HotspotWriter.Settings()
        settings.output_superstar = True
        with HotspotWriter("testdata/hs_io/minimal_all_grids", settings=settings) as w:
            w.write(a)

    def test_write_fake_multi(self):
        a = self.generate_fake(buriedness=True, superstar=True)
        b = self.generate_fake(buriedness=True, superstar=True)
        settings = HotspotWriter.Settings()
        settings.output_superstar = True
        with HotspotWriter("testdata/hs_io/minimal_multi_all_grids", settings=settings) as w:
            w.write([a, b])


class TestHotspotWriter(unittest.TestCase):
    def test_write_pymol_isosurfaces(self):
        # test out.zip prepared, generate minimal pymol commands to test isosurface gen code
        settings = HotspotWriter.Settings()
        writer = HotspotWriter("testdata/hs_io/minimal_all_grids", settings=settings)  # we won't actually write

        # pymol file initialised in the writer init function, therefore the unzip code is already in place
        writer.pymol_out.commands += writer._write_pymol_isosurfaces({"apolar": None, "donor": None, "acceptor": None},
                                                                     "hotspot",
                                                                     "hotspot",
                                                                     "fhm")

        writer.pymol_out.commands += writer._write_pymol_isosurfaces({"apolar": None, "donor": None, "acceptor": None},
                                                                     "hotspot",
                                                                     "hotspot",
                                                                     "superstar")

        writer.pymol_out.write("testdata/hs_io/minimal_all_grids/test_write_pymol_isosurfaces.py")

    def test_write_pymol_isoslider(self):
        # read in manually
        path = "testdata/hs_io/minimal_all_grids/out.zip"
        base = tempfile.mkdtemp()
        with zipfile.ZipFile(path) as hs_zip:
            hs_zip.extractall(base)

        base = os.path.join(base, "hotspot")

        interactions = ["donor", "acceptor", "apolar"]
        super_grids = {p: Grid.from_file(os.path.join(base, f"{p}.grd")) for p in interactions}
        superstar_grids = {p: Grid.from_file(os.path.join(base, f"superstar_{p}.grd")) for p in interactions}
        prot = Protein.from_file(os.path.join(base, "protein.pdb"))

        hr = Results(super_grids=super_grids,
                     protein=prot,
                     superstar=superstar_grids)

        hr.identifier = "hotspot"

        settings = HotspotWriter.Settings()
        settings.output_superstar = True

        writer = HotspotWriter("testdata/hs_io/minimal_all_grids", settings=settings)  # we won't actually write

        writer.pymol_out.commands += writer._write_pymol_isosurfaces(hr.super_grids,
                                                                     "hotspot",
                                                                     "hotspot",
                                                                     "fhm")

        writer.pymol_out.commands += writer._write_pymol_isosurfaces(hr.superstar,
                                                                     "hotspot",
                                                                     "hotspot",
                                                                     "superstar")

        writer._write_pymol_isoslider(hr)

        writer.pymol_out.write("testdata/hs_io/minimal_all_grids/test_write_pymol_isoslider.py")

    # def test_get_labels(self):
    #     labs = self.result.grid_labels()



