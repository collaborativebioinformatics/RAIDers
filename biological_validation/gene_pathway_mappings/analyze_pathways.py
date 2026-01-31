import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

from biological_validation.gene_pathway_mappings.als_master_list import (
    GENE_PATHWAY_METADATA,
    GENE_TO_PATHWAY,
    ALL_PATHWAYS
)

from biological_validation.gene_pathway_mappings.als_master_list import ALL_PATHWAYS


# Load data
df = pd.read_csv("synthetic_patients_pathways.csv")

# Map clusters to severity labels
cluster_to_severity = {
    3: 'Control',
    0: 'Mild',
    4: 'Moderate',
    1: 'Moderately Severe',
    2: 'Severe'
}

df['severity_label'] = df['severity_cluster'].map(cluster_to_severity)

severity_order = ['Control', 'Mild', 'Moderate', 'Moderately Severe', 'Severe']
df['severity_label'] = pd.Categorical(df['severity_label'], categories=severity_order)

# Basic stats
print(f'Total number of patients: {len(df)}')
print(f'\nPatients per cluster:')
print(df['severity_label'].value_counts().sort_index())

print(f'\nPatients with at least one pathway disrupted:')
pathway_cols = [c for c in df.columns if c.endswith('_binary')]
df['any_pathway'] = df[pathway_cols].sum(axis=1) > 0  # Boolean: True if > 0

print(df['any_pathway'].value_counts())

# Patients with variants but no pathway disruption
no_pathway_but_has_variants = df[(df['n_variants'] > 0) & (df['any_pathway'] == False)]
print(f"Patients with variants but no validated genes: {len(no_pathway_but_has_variants)}")

# Patients with validated genes
has_validated_genes = df[df['n_genes'] > 0]
print(f"Patients with validated genes: {len(has_validated_genes)}")

# Should match patients with pathways
print(f"Patients with pathway disruption: {df['any_pathway'].sum()}")


######## 4.1 Pathway Prevalence by Severity ####################
# Calculate pathway prevalence for each cluster
def pathway_prevalence(df):
    prevalence_results = []
    for cluster in df['severity_label'].unique():
        cluster_df = df[df['severity_label'] == cluster]
        n_patients = len(cluster_df)
        for pathway in ALL_PATHWAYS:
            pathway_patient_count = cluster_df[f'{pathway}_binary'].sum()
            pathway_prevalence = pathway_patient_count / n_patients * 100

            prevalence_results.append({
                'cluster': cluster,
                'pathway': pathway,
                'n_patients': n_patients,
                'n_affected': pathway_patient_count,
                'prevalence_pct': pathway_prevalence
            })
    return prevalence_results


# Calculate average pathway scores for each cluster
def average_pathway_scores(df):
    # Out of the affected patients, what is the average score
    average_scores_results = []
    for cluster in df['severity_label'].unique():
        cluster_df = df[df['severity_label'] == cluster]
        for pathway in ALL_PATHWAYS:
            affected_patients = cluster_df[cluster_df[f'{pathway}_binary'] == 1]
            avg_score = affected_patients[f'{pathway}_score'].mean() if len(affected_patients) > 0 else 0

            average_scores_results.append({
                'cluster': cluster,
                'pathway': pathway,
                'avg_score': avg_score
            })
    return average_scores_results

# Create heatmap: clusters (rows) * pathways (columns)
def prevalence_heatmap(heatmap_data):
    # Pivot to heatmap shape
    heatmap_matrix = heatmap_data.pivot(index='cluster', columns='pathway', values='prevalence_pct')

    # Create annotations for avg scores
    annotations = heatmap_data.pivot(index='cluster', columns='pathway', values='avg_score')

    # Create the heatmap
    plt.figure(figsize=(14, 6))  # Optional: Adjusts the figure size
    sns.heatmap(
        heatmap_matrix,
        annot=annotations,  # Annotate with Average scores
        cmap='YlOrRd',  # yellow, orange, red
        fmt=".1f",  # One decimal place for scores
        linewidths=0.5,
        linecolor='gray',
        cbar_kws={'label': 'Prevalence (%)'},  # Label the color bar
        # vmin = 0,                               # Prevalence starts at 0
        # vmax = 100                              # Prevalence max at 100%
    )

    # Add title and display the plot
    plt.title('Figure 4.1: Pathway Prevalence Across Severity Tiers', fontsize=14, fontweight='bold')
    plt.xlabel('Pathway', fontsize=12)
    plt.ylabel('Severity Cluster', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Save and show
    # plt.savefig('figure_4_1_pathway_prevalence_heatmap.png', dpi=300, bbox_inches='tight')
    # plt.show()
    plt.close()


#prevalence_df = pd.DataFrame(pathway_prevalence(df))
#avg_scores_df = pd.DataFrame(average_pathway_scores(df))
#merged_df = prevalence_df.merge(avg_scores_df, on=['cluster','pathway'], how='left')
#prevalence_heatmap(merged_df)

# Sort by cluster then pathway
#table_4_1 = merged_df.sort_values(['cluster', 'pathway']).copy()
## Round decimals for readability
#table_4_1['prevalence_pct'] = table_4_1['prevalence_pct'].round(1)
#table_4_1['avg_score'] = table_4_1['avg_score'].round(2)
# Save as CSV
#table_4_1.to_csv('table_4_1_pathway_prevalence_by_cluster.csv', index=False)
#print("✓ Table 4.1 saved to 'table_4_1_pathway_prevalence_by_cluster.csv'")
#print(f"\nTable shape: {table_4_1.shape}")
#print(f"\nFirst few rows:")
#print(table_4_1.head(10))


####  #### 4.2 Pathway Co-occurence Matrices ####################
def create_cooccurence_matrix(df):
    cooccurence_matrices = {}

    for cluster in df['severity_label'].unique():
        cluster_df = df[df['severity_label'] == cluster]


        matrix = pd.DataFrame(index=ALL_PATHWAYS, columns=ALL_PATHWAYS, dtype=float)

        for current_pathway in ALL_PATHWAYS:
            patients_with_current_pathway = cluster_df[cluster_df[f'{current_pathway}_binary'] == 1]
            n_patients_current_pathway = len(patients_with_current_pathway)


            for accompanying_pathway in ALL_PATHWAYS:
                n_patients_both_pathways = patients_with_current_pathway[f'{accompanying_pathway}_binary'].sum()

                if n_patients_current_pathway > 0:
                    cooccurence_pct = (n_patients_both_pathways / n_patients_current_pathway) * 100

                else:
                    cooccurence_pct = 0
                matrix.loc[current_pathway, accompanying_pathway] = cooccurence_pct

        # Store matrix for this cluster
        cooccurence_matrices[cluster] = matrix

    return cooccurence_matrices


def plot_cooccurence_heatmap(heatmap_matrix, cluster_name, blank_threshold = None):
    # Create display matrix (either full values or with blanking)
    display_matrix = heatmap_matrix.copy()
    if blank_threshold:
        # Blank out values below threshold for display
        display_matrix = display_matrix.map(lambda x: '' if x < blank_threshold else f'{x:.1f}')
        annot = display_matrix  # Use the blanked strings
        fmt = ''  # Don't format since we already have strings
    else:
        annot = True  # Show all values
        fmt = '.1f'  # One decimal place

    # Create the heatmap
    plt.figure(figsize=(12, 10))  # Optional: Adjusts the figure size
    sns.heatmap(
        heatmap_matrix,
        annot=annot,                      # Annotate with Average scores
        cmap='YlOrRd',                          # yellow, orange, red
        fmt=fmt,                              # One decimal place for scores
        linewidths=0.5,
        linecolor='gray',
        cbar_kws={'label': 'Co-occurence (%)'},   # Label the color bar
    )

    # Add title and display the plot
    plt.title(f'Pathway Co-occurence Matrix: {cluster_name}', fontsize=14, fontweight='bold')
    plt.xlabel('Pathway (Y)', fontsize=12)
    plt.ylabel('(Pathway (X): "Of patients with X, what % have Y?"', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save and show
    #plt.savefig('figure_4_1_pathway_prevalence_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

def extract_high_cooccurrence(matrix, cluster_name, threshold = 50):
    results = []

    for pathway_x in matrix.index:
        for pathway_y in matrix.columns:
            # Skip the diagonal
            if pathway_x == pathway_y:
                continue

            cooccurence = matrix.loc[pathway_x, pathway_y]

            if cooccurence >= threshold:
                results.append({
                    'cluster': cluster_name,
                    'pathway_X': pathway_x,
                    'pathway_Y': pathway_y,
                    'cooccurence_pct': cooccurence,
                })

    # Convert to df
    df = pd.DataFrame(results)

    # Sort df
    if not df.empty:
        df = df.sort_values(['pathway_X', 'pathway_Y', 'cooccurence_pct'], ascending = False)

    return df


pd.set_option('display.float_format', '{:.2f}'.format)
cluster_cooccurence_matrices = create_cooccurence_matrix(df)

for cluster, matrix in cluster_cooccurence_matrices.items():
    plot_cooccurence_heatmap(matrix, cluster)


all_high_cooccurrence = []

for cluster_name, matrix in cluster_cooccurence_matrices.items():
    high_pairs = extract_high_cooccurrence(matrix, cluster_name, threshold=50)
    all_high_cooccurrence.append(high_pairs)

# Combine all clusters
high_cooccurrence_df = pd.concat(all_high_cooccurrence, ignore_index=True)

print(high_cooccurrence_df[['cluster', 'pathway_X', 'pathway_Y', 'cooccurence_pct' ]])
#high_cooccurrence_df.to_csv('table_4_2_high_cooccurrence_pairs.csv', index=False)


####  #### 4.3 Compare Co-occurrence vs Correlation ####################



