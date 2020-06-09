
from os.path import join
import tempfile
import zipfile
import math
from pymol import cmd, finish_launching
from pymol.cgo import *

finish_launching()

_myarrow_cylinder = [CYLINDER, 1, 1, 0, 3.7, 4.6, 0.0, 0.1, 1, 0, 0, 0, 0, 1] 

_myarrow_cone = [CONE, 3.7, 4.6, 0.0, 4, 5, 0, 0.17500000000000002, 0, 0, 0, 1, 0, 0, 1, 1, 0.0] 

_myarrow = _myarrow_cylinder + _myarrow_cone 

cmd.load_cgo(_myarrow, "myarrow", 1)