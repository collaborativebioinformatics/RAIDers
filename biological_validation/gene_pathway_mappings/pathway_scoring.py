import pandas as pd
from biological_validation.gene_pathway_mappings.als_master_list import (
    GENE_PATHWAY_METADATA,
    GENE_TO_PATHWAY,
    ALL_PATHWAYS
)



def load_clinvar_mapping(clinvar_path, validated_genes_only=True):
    """
    Load both rsID and chr:pos mappings from Clinvar
    The 'rsID' column contains either:
    - Real rsIDs
    - chr:pos format for variants without rsIDs
    Parameters
    ----------
    clinvar_path: str
        File path to Clinvar data
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
    # Isolate rsID by removing geno_ prefix
    rsid = variant_id.replace('geno_', '')

    # Check if rsID is in mapping dictionary
    if rsid in rsid_mapping:
        # Return gene if found
        gene = rsid_mapping[rsid]
        return gene
    else:
        # Return None if rsID does  not exist in mapping
        return None


def get_pathway_scores(patient_genes, weight_primary = 1.0, weight_secondary = 0.5) -> dict:
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
    # Initialize score dictionary
    pathway_scores = {f'{pathway}_score': 0.0 for pathway in ALL_PATHWAYS}

    # for each variant gene in the patient's list
    for gene in patient_genes:
        # if gene is in metadata pathway
        if gene in GENE_PATHWAY_METADATA:
            # Get primary pathway
            prim_pathway = GENE_PATHWAY_METADATA[gene]['primary']
            pathway_scores[f'{prim_pathway}_score'] += weight_primary

            for sec_pathway in GENE_PATHWAY_METADATA[gene]['secondary']:
                pathway_scores[f'{sec_pathway}_score'] += weight_secondary

    return pathway_scores



def get_pathway_binary(pathway_scores) -> dict:
   """
   Convert pathway score to binary representation
   Parameters
   ----------
   pathway_scores: dict
   Dictionary that provides patient's scores (value) for each pathway (key) from get_pathway_scores()

   Returns
   dict:
   A dictionary that represents binary representation of pathway presence/absence
   -------
   """
   #Initialize a default dictionary for all pathways where pathway is the key and
   # value is 0 (for pathway is absent)
   pathway_binary = {f'{pathway}_binary': 0 for pathway in ALL_PATHWAYS}

   for pathway_with_score_suffix,score in pathway_scores.items():
       if score > 0:
           pathway_name = pathway_with_score_suffix.replace('_score', '')
           pathway_binary[f'{pathway_name}_binary'] = 1

   return pathway_binary