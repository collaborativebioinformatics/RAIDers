import pandas as pd
from pathway_scoring import (
    load_clinvar_mapping,
    variant_to_gene,
    get_pathway_scores,
    get_pathway_binary
)

# Load original synthetic data
df = pd.read_csv("../../nvflare/kmeans/synthetic_patients.csv")
df_severity = pd.read_csv("synthetic_patients_with_clusters.csv")
df = df.merge(df_severity[['patient_id', 'cluster']], on='patient_id', how='left')

# Get variant columns from synthetic data
variant_cols = [col for col in df.columns if col.startswith('geno')]

# Load clinvar mapping
clinvar_mapping = load_clinvar_mapping("clinvar.cleaned.csv")

results = []

# For each patient get necessary information
for _,patient_row in df.iterrows():
    # get patient_id: Unique identifier
    patient_id = patient_row["patient_id"]

    # phase1_cluster: Severity cluster (0=mild, 1, 2, 3, 4=severe) from Phase 1
    severity_cluster = patient_row["cluster"]

    #superpopulation: Ancestry (AFR, AMR, EAS, EUR, SAS)
    superpopulation = patient_row["superpopulation"]

    # Get patient variant columns
    patient_variant_cols = patient_row[variant_cols]
    # Filter for columns that have 1's
    has_variant = patient_variant_cols[patient_variant_cols == 1]
    # Get the column names
    variant_ids = has_variant.index.tolist()

    # Get variant to gene list
    patient_genes = [variant_to_gene(variant_id, clinvar_mapping) for variant_id in variant_ids]
    # Filter out None Values
    patient_genes = [gene for gene in patient_genes if gene is not None]

    # n_variants: Total number of variants this patient has
    n_variants = len(variant_ids)

    # n_unique_genes: Total unique genes
    n_unique_genes = len(set(patient_genes))

    # pathway_[name]_score: Integer (0, 1, 2, 3...) - how many genes disrupted in this pathway?
    patient_pathway_scores = get_pathway_scores(patient_genes)

    # pathway_[name]: Binary (0/1) - is this pathway disrupted?
    patient_binary_pathway = get_pathway_binary(patient_pathway_scores)

    # Store patient results in dictionary
    patient_record = {
        'patient_id': patient_id,
        'superpopulation': superpopulation,
        'severity_cluster': severity_cluster,
        'n_variants': n_variants,
        'n_genes': n_unique_genes,
    }
    patient_record.update(patient_pathway_scores)
    patient_record.update(patient_binary_pathway)

    results.append(patient_record)

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv('synthetic_patients_pathways.csv', index=False)
print(f'Saved {len(results_df)} patients to synthetic_patients_pathways.csv')