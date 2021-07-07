"""

These test will write a series of visualisation files which then need to manually be run within PyMOL
for testing purposes.

(A bit clunky but it's okay for now)

"""
import unittest
from hotspots.wrapper_pymol import PyMOLFile, PyMOLCommands


class TestPyMOLWrapper(unittest.TestCase):
    def setUp(self):
        self.pymol_file = PyMOLFile()
        self.pomegranate = (240 / 255, 52 / 255, 52 / 255, 0.5)

    def test_sphere(self):
        self.pymol_file.commands += PyMOLCommands.sphere("_mysphere",
                                                         self.pomegranate,
                                                         (0, 0, 0),
                                                         radius=2)
        self.pymol_file.commands += PyMOLCommands.load_cgo("_mysphere", "mysphere")

        self.pymol_file.commands += PyMOLCommands.sphere("_mysphere1",
                                                         self.pomegranate,
                                                         (1, 1, 1),
                                                         radius=1)
        self.pymol_file.commands += PyMOLCommands.load_cgo("_mysphere1", "mysphere1")

        self.pymol_file.write("testdata/wrapper_pymol/test_sphere.py")

    def test_cylinder(self):
        self.pymol_file.commands += PyMOLCommands.cylinder("_mycylinder",
                                                           (0, 0, 0),
                                                           (1, 1, 1),
                                                           2,
                                                           (1, 0, 0, 0.5))
        self.pymol_file.commands += PyMOLCommands.load_cgo("_mycylinder", "mycylinder")

        self.pymol_file.commands += PyMOLCommands.cylinder("_mycylinder2",
                                                           (0, 0, 0),
                                                           (10, 10, 10),
                                                           0.2,
                                                           (1, 0, 0, 1),
                                                           (0.5, 0.5, 0))
        self.pymol_file.commands += PyMOLCommands.load_cgo("_mycylinder2", "mycylinder2")

        self.pymol_file.write("testdata/wrapper_pymol/test_cylinder.py")

    def test_cone(self):
        self.pymol_file.commands += PyMOLCommands.cone("_mycone",
                                                       (0, 0, 0),
                                                       (10, 10, 10),
                                                       5,
                                                       (1, 0, 0))
        self.pymol_file.commands += PyMOLCommands.load_cgo("_mycone", "mycone")

        self.pymol_file.write("testdata/wrapper_pymol/test_cone.py")

    def test_arrow(self):
        self.pymol_file.commands += PyMOLCommands.arrow("_myarrow",
                                                        (1, 1, 0),
                                                        (4, 5, 0),
                                                        (1, 0, 0),
                                                        (0, 0, 1))
        self.pymol_file.commands += PyMOLCommands.load_cgo("_myarrow", "myarrow")

        self.pymol_file.write("testdata/wrapper_pymol/test_arrow.py")

    def test_pseudoatom(self):
        self.pymol_file.commands += PyMOLCommands.pseudoatom("mypseudoatom",
                                                             (1, 1, 1),
                                                             self.pomegranate)

        self.pymol_file.commands += PyMOLCommands.pseudoatom("mypseudoatom2",
                                                             (2, 2, 2),
                                                             self.pomegranate,
                                                             "31.2")

        self.pymol_file.write("testdata/wrapper_pymol/test_pseudoatom.py")

    def test_line(self):
        self.pymol_file.commands += PyMOLCommands.line("myline",
                                                       (1, 1, 1),
                                                       (2, 2, 2),
                                                       self.pomegranate,
                                                       width=4)

        self.pymol_file.commands += PyMOLCommands.line("myline2",
                                                       (1, 1, 1),
                                                       (0, 0, 0),
                                                       self.pomegranate,
                                                       width=10)

        self.pymol_file.write("testdata/wrapper_pymol/test_line.py")

    def test_grid(self):
        self.pymol_file.commands += PyMOLCommands.load("test_grid.grd", "mygrid")
        self.pymol_file.write("testdata/wrapper_pymol/test_grid.py")

    def test_isosurface(self):
        self.pymol_file.commands += PyMOLCommands.set_color("pomegranate", self.pomegranate)
        self.pymol_file.commands += PyMOLCommands.isosurface("test_grid.grd",
                                                             "mygrid",
                                                             "mysurface",
                                                             4,
                                                             "pomegranate")

        self.pymol_file.write("testdata/wrapper_pymol/test_isosurface.py")

    def test_load(self):
        self.pymol_file.commands += PyMOLCommands.load("test_protein.pdb", "myprotein")
        self.pymol_file.write("testdata/wrapper_pymol/test_load.py")

    def test_group(self):
        self.pymol_file.commands += PyMOLCommands.line("myline1",
                                                       (0, 0, 0),
                                                       (1, 1, 1),
                                                       self.pomegranate,
                                                       width=10)

        self.pymol_file.commands += PyMOLCommands.pseudoatom("mypseudoatom1",
                                                             (0, 0, 0),
                                                             self.pomegranate)

        self.pymol_file.commands += PyMOLCommands.pseudoatom("mypseudoatom2",
                                                             (1, 1, 1),
                                                             self.pomegranate)
        self.pymol_file.commands += PyMOLCommands.group("mygroup",
                                                        ["myline1", "mypseudoatom1", "mypseudoatom2"])

        self.pymol_file.write("testdata/wrapper_pymol/test_group.py")




