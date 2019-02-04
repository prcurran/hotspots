# Welcome to the fragment-hotspot-maps wiki!


![fragment hotspots](http://fragment-hotspot-maps.ccdc.cam.ac.uk/static/cover_small.jpg)


***************
Hotspots API
***************

The Hotspots API is the python package for the Fragment Hotspot Maps project, a knowledge-based method for determining small molecule binding "hotspots" 

For more information on this method: 
    -  - Radoux, C.J. et. al., Identifying the Interactions that Determine Fragment Binding at Protein Hotspots J. Med. Chem. 2016, 59 (9), 4314-4325 [dx.doi.org/10.1021/acs.jmedchem.5b01980]


***************
Getting Started
***************

Although the Hotspots API is publicly available under the MIT license, it is dependant on the CSD python API - a commercial package. 
If you are an academic user, it's likely your institution will have a license. If you are unsure if you have a 
license or would like to enquire about purchasing one, please contact support@ccdc.cam.ac.uk


Please note, this is an academic project and we would therefore welcome feedback, contributions and collaborations.
If you have any queries regarding this package please contact us (pcurran@ccdc.cam.ac.uk)!


NB: We recommend installing on a Linux machine

===================
Installation
===================

-------------------------
Step 1: Install CSDS 2019
-------------------------
Available from CCDC downloads page `here <https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/>`_.


You will need a valid site number and confirmation code, this will have been emailed to you when you bought your CSDS 2019 license

You may need to set the following environment variables:

.. code-block:: shell

    export CSDBASE="<path_to_CSDS_installation>"
    export CSDHOME=%CSDBASE%/CSD_2019
    export SUPERSTAR_ROOT=%CSDBASE%/Discovery_2019/SuperStar/
    export SUPERSTAR_ISODIR=%CSDBASE%/Discovery_2019/GOLD
    export GOLD_DIR=%CSDHOME%/isostar_files/istr/



-------------------------
Step 2: Install Ghecom
-------------------------
Available from Ghecom download page `here <http://strcomp.protein.osaka-u.ac.jp/ghecom/download_src.html>`_.

"The source code of the ghecom is written in C, and developed and executed on
the linux environment (actually on the Fedora Core).  For the installation,
you need the gcc compiler.  If you do not want to use it, please change the
"Makefile" in the "src" directory."

Download the file "ghecom-src-[date].tar.gz" file.

.. code-block:: shell

    tar zxvf ghecom-src-[date].tar.gz
    cd src
    make
    export GHECOM_EXE="<download_directory>"
	
	
------------------------------------------------
Step 3: Create a conda environment (recommended)
------------------------------------------------

.. code-block:: shell

    conda create -n hotspots_env python=2.7

------------------------------------------------
Step 4: Install RDKit and CSD python API
------------------------------------------------
Download the standalone CSD python API package from `here <https://www.ccdc.cam.ac.uk/forum/csd_python_api/General/06004d0d-0bec-e811-a889-005056977c87>`_.

.. code-block:: shell

	conda install -c rdkit -n hotspots_env rdkit
	conda install csd-python-api-2.x.x-linux-py2.7-conda.tar.bz2
	
------------------------------------------------
Step 5: Install Hotspots API
------------------------------------------------

.. code-block:: shell

    source activate hotspots_env
    pip install hotspots


... and you're ready to go!

*********************
Running a Calculation
*********************

===================
Protein Preparation
===================

The first step is to make sure your protein is correctly prepared for the calculation. The structures should be
protonated with small molecules and waters removed. Any waters or small molecules left in the structure will be included
in the calculation.

One way to do this is to use the CSD Python API:


.. code-block:: python
    
    from ccdc.protein import Protein

    prot = Protein.from_file('protein.pdb')
    prot.remove_all_waters()
    prot.add_hydrogens()
    for l in prot.ligands:
        prot.remove_ligand(l.identifier)


For best results, manually check proteins before submitting them for calculation.


=================================
Calculating Fragment Hotspot Maps
=================================

Once the protein is prepared, the :class:`hotspots.calculation.Runner` object can be used to perform the calculation:

.. code-block:: python

    from hotspots.calculation import Runner

    r = Runner()
    results = Runner.from_protein(prot)
	

Alternatively, for a quick calculation, you can supply a PDB code and we will prepare the protein as described above:

.. code-block:: python

    r = Runner()
    results = Runner.from_pdb("1hcl")

============================
Reading and Writing Hotspots
============================

-------
Writing
-------

The :mod:`hotspots.hs_io` module handles the reading and writing of both :class:`hotspots.calculation.results`
and :class:`hotspots.best_volume.Extractor` objects. The output `.grd` files can become quite large, but are highly
compressible, therefore the results are written to a `.zip` archive by default, along with a PyMOL run script to
visualise the output.

.. code-block:: python

    from hotspots import hs_io
	
    out_dir = "results/pdb1"

    # Creates "results/pdb1/out.zip"
    with HotspotWriter(out_dir, grid_extension=".grd", zip_results=True) as w:
        w.write(results)

-------
Reading
-------

If you want to revisit the results of a previous calculation, you can load the `out.zip` archive directly into a
:class:`hotspots.calculation.results` instance:

.. code-block:: python

    from hotspots import hs_io

    results = hs_io.HotspotReader('results/pdb1/out.zip').read()

****************
Using the Output
****************

While Fragment Hotspot Maps provide a useful visual guide, the grid-based data can be used in other SBDD analysis.

=======
Scoring
=======
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
    


To learn about other ways you can use the Hotspots API please read our API documentation.
