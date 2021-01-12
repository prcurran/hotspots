# Unfortunately RDKit and matplotlib are both linked to specific underlying 
# versions of compiled libraries. The CCDC code is less fussy, but this means if
# these modules are imported after CCDC modules the imports of rdkit & matplolib 
# can fail depending on the exact version of the underlying library pulled in
# if it happens the version the CCDC import pulled in was one that they dont work with.
#
from rdkit import *
from rdkit.Chem.rdmolops import *
from matplotlib import *

from hotspots import calculation


__author__ = "Chris Radoux, Peter Curran, Mihaela Smilova"
__copyright__ = None
__credits__ = None
__license__ = None
__version__ = "1.0.3"
__maintainer__ = "Peter Curran"
__email__ = "pcurran@ccdc.cam.ac.uk"
__status__ = "Development"
name = "hotspots"
