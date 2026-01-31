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


####  #### 4.2A Pathway Co-occurence Matrices ####################
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
    #plt.savefig(f'figure_4_1_pathway_cooccurrence_heatmap_{cluster_name}.png', dpi=300, bbox_inches='tight')
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
cooccurence_matrices = create_cooccurence_matrix(df)

#for cluster, matrix in cluster_cooccurence_matrices.items():
#    print(f'Cluster {cluster}: {matrix}')
#for cluster, matrix in cluster_cooccurence_matrices.items():
    #plot_cooccurence_heatmap(matrix, cluster)


#all_high_cooccurrence = []

#for cluster_name, matrix in cluster_cooccurence_matrices.items():
#    high_pairs = extract_high_cooccurrence(matrix, cluster_name, threshold=50)
#    all_high_cooccurrence.append(high_pairs)

# Combine all clusters
#high_cooccurrence_df = pd.concat(all_high_cooccurrence, ignore_index=True)

#print(high_cooccurrence_df[['cluster', 'pathway_X', 'pathway_Y', 'cooccurence_pct' ]])
#high_cooccurrence_df.to_csv('table_4_2_high_cooccurrence_pairs.csv', index=False)


######## 4.2B Pathway Correlation Matrices ####################

def create_correlation_matrix(df, variance_threshold = 0.01):
    correlation_matrices = {}
    for cluster in df['severity_label'].unique():
        cluster_df = df[df['severity_label'] == cluster]

        pathways_with_variance = []
        for pathway in ALL_PATHWAYS:
            score_col = f'{pathway}_score'
            variance = cluster_df[score_col].var()

            if variance > variance_threshold:
                pathways_with_variance.append(pathway)

        # Step 2: Calculate correlation only for those pathways
        if len(pathways_with_variance) > 1:  # Need at least 2 pathways to correlate
            score_columns = [f'{p}_score' for p in pathways_with_variance]
            pathway_scores = cluster_df[score_columns]

            matrix = pathway_scores.corr(method='pearson')

            matrix.index = pathways_with_variance
            matrix.columns = pathways_with_variance
        else:
            matrix = pd.DataFrame() # Empty

        correlation_matrices[cluster] = matrix
    return correlation_matrices

def plot_correlation_heatmap(correlation_matrix, cluster_name, blank_threshold = None):
    if correlation_matrix.empty or len(correlation_matrix) < 2:
        print(f'Skipping {cluster_name}: insufficient data for correlation heatmap')
        return
    # Create display matrix (either full values or with blanking)
    display_matrix = correlation_matrix.copy()
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
        correlation_matrix,
        annot=annot,  # Annotate with Average scores
        cmap='RdBu_r',  # Red-Blue Reversed (red = pos, blue = negative)
        fmt=fmt,  # One decimal place for scores
        linewidths=0.5,
        linecolor='gray',
        cbar_kws={'label': 'Correlation'},  # Label the color bar
        vmin = -1,
        vmax = 1
    )

    # Add title and display the plot
    plt.title(f'Pathway Correlation Matrix: {cluster_name}', fontsize=14, fontweight='bold')
    plt.xlabel('Pathway (X)', fontsize=12)
    plt.ylabel('Pathway (Y)', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save and show
    plt.savefig(f'figure_4_3_pathway_correlation_heatmap_{cluster_name}.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()

correlation_matrices = create_correlation_matrix(df)
#for cluster, matrix in correlation_matrices.items():
#    print(f'\n{cluster}: {matrix}')

#for cluster, matrix in correlation_matrices.items():
#    plot_correlation_heatmap(matrix, cluster)

######## 4.3 Compare Co-occurrence vs Correlation ####################
def combine_cooccurrence_correlation(cooccurrence_matrices, correlation_matrices):
    """
    Combine cooccurrence and correlation data for scatter plot analysis
    Parameters
    ----------
    cooccurrence_matrices
    correlation_matrices

    Returns
    Dataframe
    -------

    """
    results = []

    for cluster in cooccurrence_matrices.keys():
        cooccur_matrix = cooccurrence_matrices[cluster]
        corr_matrix = correlation_matrices[cluster]

        if cooccur_matrix.empty or corr_matrix.empty:
            continue

        common_pathways = list(set(cooccur_matrix.index) & set(corr_matrix.index))

        for pathway_X in common_pathways:
            for pathway_Y in common_pathways:
                if pathway_X == pathway_Y:
                    continue

                cooccur_value = cooccur_matrix.loc[pathway_X, pathway_Y]
                corr_value = corr_matrix.loc[pathway_X, pathway_Y]

                results.append({
                    'cluster': cluster,
                    'pathway_X': pathway_X,
                    'pathway_Y': pathway_Y,
                    'cooccurrence_pct': cooccur_value,
                    'correlation': corr_value
                })

    return pd.DataFrame(results)

def plot_cooccurrence_vs_correlation(combined_df, cluster_name):
    cluster_data = combined_df[combined_df['cluster'] == cluster_name]

    # Skip if no data
    if cluster_data.empty:
        print(f'Skipping {cluster_name}: no overlapping pathway pairs')
        return
    plt.figure(figsize=(10, 8))

    # Color-code points by pattern type
    colors = []
    for _, row in cluster_data.iterrows():
        cooccur = row['cooccurrence_pct']
        corr = row['correlation']

        if cooccur > 50 and corr > 0.5:
            colors.append('green')  # Upper right: Dose-dependent
        elif cooccur > 50 and abs(corr) < 0.3:
            colors.append('orange')  # Lower right: Threshold effect
        elif cooccur < 30 and corr < -0.3:
            colors.append('red')  # Lower left: Distinct subtypes
        else:
            colors.append('lightgray')  # Other


    plt.scatter(
        x=cluster_data['cooccurrence_pct'],
        y=cluster_data['correlation'],
        c = colors,
        alpha=0.6,
        s=80,
    )

    plt.xlabel('Co-occurence (%)', fontsize=12)
    plt.ylabel(' Pearson Correlation', fontsize=12)
    plt.title(f'Co-occurrence vs Correlation: {cluster_name}', fontsize=14, fontweight='bold')

    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    plt.axvline(x=50, color='gray', linestyle='--', linewidth=0.5)

    # Add quadrant labels
    plt.text(75, 0.85, 'Dose-dependent\n(cascading failure)',
             fontsize=10, ha='center', style='italic', color='darkgreen', weight='bold')
    plt.text(75, -0.5, 'Threshold effect\n(binary disruption)',
             fontsize=10, ha='center', style='italic', color='darkorange', weight='bold')
    plt.text(25, 0.85, 'Rare patterns',
             fontsize=9, ha='center', style='italic', color='gray')
    plt.text(25, -0.5, 'Independent/\nDistinct subtypes',
             fontsize=9, ha='center', style='italic', color='gray')

    plt.xlim(0,100)
    plt.ylim(-1,1)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='green', edgecolor='white', label='Dose-dependent (co-occ >50%, r >0.5)'),
        Patch(facecolor='orange', edgecolor='white', label='Threshold effect (co-occ >50%, |r| <0.3)'),
        Patch(facecolor='red', edgecolor='white', label='Distinct subtypes (co-occ <30%, r <-0.3)'),
        Patch(facecolor='lightgray', edgecolor='white', label='Other patterns')
    ]
    plt.legend(handles=legend_elements, loc='lower left', fontsize=10, framealpha=0.9)

    # Grid
    plt.grid(True, alpha=0.2, linestyle=':')

    plt.tight_layout()
    plt.savefig(f'figure_4_4_cooccurence_vs_correlation_{cluster_name}.png', dpi=300)
    plt.show()
    plt.close()

combined_matrices = combine_cooccurrence_correlation( cooccurence_matrices ,correlation_matrices)
print(combined_matrices)

for cluster_name in combined_matrices['cluster'].unique():
    plot_cooccurrence_vs_correlation(combined_matrices, cluster_name)