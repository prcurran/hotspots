************
# Hotspots API
************

[![Documentation Status](https://readthedocs.org/projects/hotspots/badge/?version=latest)](https://hotspots.readthedocs.io/en/latest/?badge=latest)	
[![License](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://github.com/prcurran/fragment_hotspot_maps/blob/master/LICENSE)	
[![Total alerts](https://img.shields.io/lgtm/alerts/g/prcurran/fragment_hotspot_maps.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/prcurran/fragment_hotspot_maps/alerts/)	
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/prcurran/fragment_hotspot_maps.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/prcurran/fragment_hotspot_maps/context:python)	
[![Language grade: JavaScript](https://img.shields.io/lgtm/grade/javascript/g/prcurran/fragment_hotspot_maps.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/prcurran/fragment_hotspot_maps/context:javascript)
[![Gitter chat](https://badges.gitter.im/gitterHQ/gitter.png)](https://gitter.im/fragmenthotspots/community)


 ![fragment hotspots](http://fragment-hotspot-maps.ccdc.cam.ac.uk/static/cover_small.jpg)


The Hotspots API is the Python package for the Fragment Hotspot Maps project,
a knowledge-based method for determining small molecule binding "hotspots".

For more information on this method:

    Radoux, C.J. et. al., Identifying the Interactions that Determine Fragment
    Binding at Protein Hotspots J. Med. Chem. 2016, 59 (9), 4314-4325
    [dx.doi.org/10.1021/acs.jmedchem.5b01980]


Getting Started
===============

Although the Hotspots API is publicly available, it is dependant on the CSD
Python API - a commercial package.

If you are an academic user, it's likely your institution will have a license.
If you are unsure if you have a license or would like to enquire about
purchasing one, please contact support@ccdc.cam.ac.uk.

Please note, this is an academic project and we would therefore welcome
feedback, contributions and collaborations. If you have any queries regarding
this package please contact us (pcurran@ccdc.cam.ac.uk)!


Installation
============


*1* Install CSDS 2019
----------------------

Available from `CCDC downloads page <https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/>`_.

You will need a valid site number and confirmation code, this will have been
emailed to you when you bought your CSDS 2019 license.


*2* Install GHECOM
-------------------

Available from `GHECOM download page <http://strcomp.protein.osaka-u.ac.jp/ghecom/download_src.html>`_.

    The source code of the GHECOM is written in C, and developed and executed on
    the linux environment (actually on the Fedora Core).  For the installation,
    you need the gcc compiler.  If you do not want to use it, please change the
    "Makefile" in the "src" directory.

Download the file ``ghecom-src-[date].tar.gz`` file.

.. code-block:: shell
    
    tar zxvf ghecom-src-[date].tar.gz
    cd src
    make
    # The executable will be located at the parent directory ..


*3* Setup an Anaconda environment (recommended)
------------------------------------------------

.. code-block:: shell
    
    # create environment and install requirements
    conda create -n hotspots python=2.7
    conda install -n hotspots \
        numpy==1.15.4 matplotlib==2.2.3 scikit-image==0.14.2 \
        pandas==0.24.1 futures==3.2.0 cython==0.29.5 tqdm==4.31.1 \
        xmltodict==0.12.0 \
        scipy scikit-learn
    
    # install the Python CSD API
    conda install -n hotspots qt==5.9.7 rdkit==2018.09.1
    conda install -n hotspots csd-python-api-2.x.x-linux-py2.7-conda.tar.bz2
    
    # install Hotspots v1.0.0
    conda run -n hotspots pip install https://github.com/prcurran/hotspots/archive/v1.0.0.zip


## Hotspots API Usage
---------------------

Start activating your Anaconda environment and setting some variables.

.. code-block:: shell

    conda activate hotspots
    export GHECOM_EXE=<path_to_GHECOM_executable>
    export CSDHOME=<path_to_CSDS_installation>/CSD_2019


## Running a Calculation
---------------------

### Protein Preparation

The first step is to make sure your protein is correctly prepared for the
calculation. The structures should be protonated with small molecules and
waters removed. Any waters or small molecules left in the structure will
be included in the calculation.

One way to do this is to use the CSD Python API:


.. code-block:: python
    
    from ccdc.protein import Protein

    prot = Protein.from_file('protein.pdb')
    prot.remove_all_waters()
    prot.add_hydrogens()
    for l in prot.ligands:
        prot.remove_ligand(l.identifier)


For best results, manually check proteins before submitting them for calculation.


### Calculating Fragment Hotspot Maps
---------------------


Once the protein is prepared, the `hotspots.calculation.Runner` object can be
used to perform the calculation:

.. code-block:: python

    from hotspots.calculation import Runner

    runner = Runner()
    results = runner.from_pdb(prot, nprocesses=11)
	

Alternatively, for a quick calculation, you can supply a PDB code and we will
prepare the protein as described above:

.. code-block:: python

    runner = Runner()
    results = runner.from_pdb("1hcl", nprocesses=11)


## Reading and Writing Hotspots
----------------------------

### Writing

The  `hotspots.hs_io` module handles the reading and writing of both  `hotspots.calculation.results
and  `hotspots.best_volume.Extractor` objects. The output `.grd` files can become quite large,
but are highly compressible, therefore the results are written to a `.zip` archive by default,
along with a PyMOL run script to visualise the output.

.. code-block:: python

    from hotspots.hs_io import HotspotWriter
	
    out_dir = "results/pdb1"

    # Creates "results/pdb1/out.zip"
    with HotspotWriter(out_dir) as writer:
        writer.write(results)

### Reading


If you want to revisit the results of a previous calculation, you can load the
`out.zip` archive directly into a `hotspots.calculation.results` instance:

.. code-block:: python

    from hotspots.hs_io import HotspotReader

    results = HotspotReader('results/pdb1/out.zip').read()



## Using the Output
---------------------

While Fragment Hotspot Maps provide a useful visual guide, the grid-based data
can be used in other SBDD analysis.

### Scoring
---------------------

One example is scoring atoms of either proteins or small molecules.

This can be done as follows: 

.. code-block:: python

    from ccdc.protein import Protein
    from ccdc.io import MoleculeReader, MoleculeWriter
    from hotspots.calculation import Runner
	
	r = Runner()
	prot = Protein.from_file("1hcl.pdb")    # prepared protein
	hs = r.from_protein(prot)
	
	# score molecule
	mol = MoleculeReader("mol.mol2")
	scored_mol = hs.score(mol)
	with MoleculeWriter("score_mol.mol2") as w:
	    w.write(scored_mol)
		
	# score protein
	scored_prot = hs.score(hs.prot)
	with MoleculeWriter("scored_prot.mol2") as w:
	    w.write(scored_prot)
    

To learn about other ways you can use the Hotspots API please see the examples
directory and read our API documentation.
