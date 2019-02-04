"""
Package data
"""
import os

def prompts():
    """
    Scrapped from RCSB PDB website- use for creating query prompts
    :return: dict
    """
    return {"HoldingsQuery": "All/Experimental Type/Molecule Type", "StructureIdQuery": "PDB ID(s)",
            "EntityIdQuery": "Entity ID(s)", "ChainIdQuery": "Chain ID(s)", "PubmedIdQuery": "PubMed ID(s)",
            "UpAccessionIdQuery": "UniProtKB Accession Number(s)", "AdvancedKeywordQuery": "Text Search",
            "TokenKeywordQuery": "mmCIF Keyword Search (Classification)", "PfamIdQuery": "Pfam Accession Number(s)",
            "UniprotGeneNameQuery": "UniProt Gene Name", "SequenceClusterQuery": "Sequence Cluster Name",
            "StructTitleQuery": "Structure Title", "StructDescQuery": "Structure Description",
            "MoleculeNameQuery": "Macromolecule Name", "PathwayQuery": "Pathway Reaction Name",
            "LargeStructureQuery": "Large Structures", "AdvancedAuthorQuery": "Author Name",
            "DepositDateQuery": "Deposit Date", "ReleaseDateQuery": "Release Date", "ReviseDateQuery": "Revision Date",
            "LastLoadQuery": "Latest Released Structures", "ModifiedStructuresQuery": "Latest Modified Structures",
            "SGProjectQuery": "Structural Genomics Project", "ChainTypeQuery": "Macromolecule Type",
            "NumberOfChainsQuery": "Number of Chains (Asymmetric Unit)",
            "BiolUnitQuery": "Number of Chains (Biological Assembly)", "NumberOfEntitiesQuery": "Number of Entities",
            "StoichiometryQuery": "Protein Stoichiometry", "PointGroupQuery": "Protein Symmetry",
            "TreeQuery_20": "Protein Symmetry Browser", "ModelCountQuery": "Number of Models",
            "CloseContactsQuery": "Number of Disulfide Bonds", "LinkConnectionQuery": "Link records",
            "MolecularWeightQuery": "Molecular Weight (Structure)",
            "SecondaryStructureQuery": "Secondary Structure Content",
            "SecondaryStructureLengthQuery": "Secondary Structure Length",
            "TreeQuery_11": "SCOP Classification Browser",
            "TreeQuery_12": "CATH Classification Browser", "TreeQueryExpression_1": "Taxonomy Browser",
            "SequenceQuery": "Sequence (BLAST/PSI-BLAST)", "WildTypeProteinQuery": "Wild Type Protein",
            "MutationQuery": "Mutation", "BlastXQuery": "Translated Nucleotide Sequence (BLASTX)",
            "MotifQuery": "Sequence Motif", "SequenceLengthQuery": "Chain Length",
            "ProteinModificationsQuery": "Protein Modifications", "TreeQuery_17": "Protein Modification Browser",
            "TreeEntityQuery_13": "Genome Location Browser", "ChemCompNameQuery": "Chemical Name",
            "ChemCompIdQuery": "Chemical ID(s)", "ChemCompDescriptorQuery": "InChI Descriptor",
            "ChemSmilesQuery": "Chemical structure (SMILES)",
            "ChemFormulaWeightQuery": "Molecular Weight (Chemical component)",
            "ChemCompFormulaQuery": "Chemical Formula", "ChemCompTypeQuery": "Chemical Component Type",
            "BindingAffinityQuery": "Binding Affinity", "NoLigandQuery": "Has Ligand(s)",
            "NoModResQuery": "Has Modified Residue(s)", "ChemCompSubCompQuery": "Sub-components",
            "BirdQuery": "Biologically Interesting Molecules (from BIRD)",
            "TreeEntityQuery_1": "Source Organism Browser (NCBI)", "ExpressionOrganismQuery": "Expression Organism",
            "TreeEntityQuery_3": "Enzyme Classification Browser", "EnzymeClassificationQuery": "Enzyme Classification",
            "TreeEntityQuery_4": "Biological Process Browser (GO)", "TreeEntityQuery_5": "Cell Component Browser (GO)",
            "TreeEntityQuery_6": "Molecular Function Browser (GO)",
            "TreeQuery_16": "Transporter Classification Browser",
            "ExpTypeQuery": "Experimental Method", "ResolutionQuery": "X-ray Resolution",
            "AverageBFactorQuery": "X-ray Average B Factor", "XrayRefinementQuery": "Refinement R Factors",
            "XrayDiffrnSourceQuery": "Diffraction Source",
            "MethodToDetermineStructQuery": "Structure Determination Method", "XrayReflnsQuery": "Reflections",
            "XrayCellQuery": "Cell Dimensions", "SoftwareQuery": "Software", "SpaceGroupQuery": "Space Group",
            "CrystalQuery": "Crystal Properties", "EmAssemblyQuery": "EM Assembly", "DetectorTypeQuery": "Detector",
            "CitationTwoQuery": "Citation", "TreeQuery_15": "Medical Subject Headings Browser",
            "PubMedQuery": "PubMed Abstract", "ExternalLinkQuery": "Has External Links"}


def common_solvents():
    """
    common solvents/ other junk in PDB entries
    :return:
    """
    f = open(os.path.join(os.path.split(__file__)[0], "excluded_het_id.txt"), 'r')
    return f.read().splitlines()
