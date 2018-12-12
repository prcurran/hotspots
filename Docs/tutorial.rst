########
Tutorial
########

This section will introduce the main functionality of the Hotspots API

***************
Getting Started
***************

*********************
Running a Calculation
*********************

===================
Protein Preparation
===================

The first step is to make sure your protein is correctly prepared for the calculation. The structures should be
protonated with small molecules and waters removed. Any waters or small molecules left in the structure will be included
in the calculation.

This can be done with the CSD Python API::

    from ccdc import Protein

    prot = Protein.from_file('protein.pdb')
    prot.remove_all_waters()
    prot.add_hydrogens()
    for l in prot.ligands:
         prot.remove_ligand(l.identifier)

=================================
Calculating Fragment Hotspot Maps
=================================

Once the protein is prepared, the :class:`hotspots.calculation.Runner` object can be used to perform the calculation::

    from hotspots.calculation import Runner

    r = Runner()
    results = Runner.from_protein(prot)

============================
Reading and Writing Hotspots
============================

-------
Writing
-------

The :mod:`hotspots.hs_io` module handles the reading and writing of both :class:`hotspots.calculation.results`
and :class:`hotsptots.best_volume.Extractor` objects. The output `.grd` files can become quite large, but are highly
compressible, therefore the results are written to a `.zip` archive by default, along with a PyMOL run script to
visualise the output. ::

    from hotspots import hs_io

    out_dir = "results/pdb1"

    # Creates "results/pdb1/out.zip"
    w = hs_io.HotspotWriter(out_dir, grid_extension=".grd", zip_results=True)
    w.write(results)


-------
Reading
-------

If you want to revisit the results of a previous calculation, you can load the `out.zip` archive directly into a
:class:`hotspots.calculation.results` instance::

    from hotspots import hs_io

    results = hs_io.HotspotReader('results/pdb1/out.zip').read()

****************
Using the Output
****************

While Fragment Hotspot Maps provide a useful visual guide, the grid-based data can be used

=======
Scoring
=======


---------------
Small molecules
---------------

--------
Proteins
--------


=======================
Creating Pharmacophores
=======================



=========================================
Extracting Hotspots and Drug-like Volumes
=========================================

------
Global
------

----
Seed
----
