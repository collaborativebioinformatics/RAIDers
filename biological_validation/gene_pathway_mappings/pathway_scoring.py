import pandas as pd
from biological_validation.gene_pathway_mappings.als_master_list import (
    GENE_PATHWAY_METADATA,
    GENE_TO_PATHWAY,
    ALL_PATHWAYS
)

def load_clinvar_mapping(clinvar_path, validated_genes_only=True):
    """
    Load clinvar csv file and convert rsID and genes into a dictionary
    Parameters
    ----------
    clinvar_path: str
    Path to clinvar csv file
    validated_genes_only: bool, default = True
    if True, only include genes from literature-validated list

    Returns:
        dict: {rsID: gene_symbol}
    -------
    """
    # Read the clinvar csv file
    df = pd.read_csv(clinvar_path)

    # Create initial mapping
    rsid_to_gene = dict(zip(df['rsID'], df['gene']))

    # Get literature validated genes
    if validated_genes_only:
        validated_genes = set(GENE_TO_PATHWAY.keys())

        # Filter clinvar csv files to only include lit validated genes
        rsid_to_gene_filtered = {
            rsid: gene
            for rsid, gene in rsid_to_gene.items()
            if gene in validated_genes
        }

        removed_count = len(rsid_to_gene) - len(rsid_to_gene_filtered)

        print(f"Total rsIDs in Clinvar: {len(rsid_to_gene)}")
        print(f"Validated genes: {len(validated_genes)}")
        print(f"rsIDs mapping to validated genes: {len(rsid_to_gene_filtered)}")
        print(f"Filtered out {removed_count} rsIDs mapping to non-validated genes")

        return rsid_to_gene_filtered
    else:
        print(f"Loaded {len(rsid_to_gene)} rsIDs -> gene mappings (unfiltered)")
        return rsid_to_gene

def variant_to_gene(variant_id, rsid_mapping):
    """
    Convert variant ID to gene symbol
    Parameters
    ----------
    variant_id
    rsid_mapping

    Returns
    -------

    """
    pass
def get_pathway_scores(patient_genes, weight_primary = 1.0, weight_secondary = 0.5)
    """
    Calculate weighted pathway score
    Parameters
    ----------
    patient_genes
    weight_primary
    weight_secondary

    Returns
    -------

    """
    pass
def get_pathway_binary(pathway_scores):
    """
    Convert scores to binary
    Parameters
    ----------
    pathway_scores

    Returns
    -------

    """
#load_clinvar_mapping("clinvar.cleaned.csv")