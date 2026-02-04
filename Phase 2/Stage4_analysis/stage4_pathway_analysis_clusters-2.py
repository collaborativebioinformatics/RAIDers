#!/usr/bin/env python3
"""
Stage 4: Pathway Co-occurrence and Correlation Analysis (CLUSTER VERSION)
==========================================================================
Analyzes pathway patterns within Phase 1 clusters (0, 1, 2, 4).

REQUIRES: Input CSV must have 'cluster' column already assigned.
          Use patients_with_pathways_weighted_clusters.csv

FIXES IN THIS VERSION:
  - Added effect sizes (Cramér's V, epsilon squared, odds ratios)
  - Fixed p-value underflow by reporting -log10(p) scale
  - Added Fisher's exact test with odds ratios for co-occurrence
  - Added confidence intervals for effect sizes
  - Better handling of extreme p-values

Analyses:
  4.1: Pathway prevalence by severity
  4.2A: Co-occurrence matrices (binary) + statistical significance
  4.2B: Correlation matrices (scores)
  4.3: Frequency vs intensity comparison
  4.4: Network visualization data (with effect sizes)
  4.5: Cross-severity comparison (with effect sizes)
"""

import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import json
import os
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Configuration
INPUT_FILE = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/ALS_Synthetic_Data/phase2_output_weighted_clusters/patients_with_pathways_weighted_clusters.csv")
OUTPUT_DIR = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/Stage4_analysis")


PATHWAYS = [
    'Proteostasis',
    'RNA_Metabolism', 
    'Cytoskeletal_Axonal_Transport',
    'Mitochondrial',
    'Excitotoxicity',
    'Vesicle_Trafficking',
    'DNA_Damage',
]

SEVERITY_ORDER = ['Unaffected', 'Mild', 'Moderate', 'Severe']

# =============================================================================
# PHASE 1 CLUSTER DEFINITIONS
# =============================================================================
# Clusters derived from Federated K-means in Phase 1
# Cluster assignments are deterministic based on severity_count ranges

CLUSTER_DEFINITIONS = {
    0: (4.0, 5.7, 'Cluster_0_Mild', 'Mild - DNA Damage pathway dominant'),
    1: (6.9, 7.8, 'Cluster_1_Moderate_SOD1', 'Moderate - SOD1/Proteostasis dominant'),
    2: (8.4, 10.0, 'Cluster_2_Severe', 'Severe - Multi-pathway involvement'),
    3: (0.0, 0.0, 'Cluster_3_Unaffected', 'Non-carriers'),
    4: (6.0, 6.6, 'Cluster_4_Moderate_ALS2', 'Moderate - ALS2/Vesicle dominant'),
}

CLUSTER_LABELS = {
    0: 'Cluster_0_Mild',
    1: 'Cluster_1_Moderate_SOD1', 
    2: 'Cluster_2_Severe',
    3: 'Cluster_3_Unaffected',
    4: 'Cluster_4_Moderate_ALS2',
}

# Carrier clusters only (excluding unaffected)
CARRIER_CLUSTERS = [0, 1, 2, 4]



# =============================================================================
# EFFECT SIZE AND P-VALUE UTILITY FUNCTIONS
# =============================================================================

def safe_log10_pvalue(p):
    """
    Convert p-value to -log10 scale, handling underflow.
    
    Returns:
        -log10(p), capped at 300 for p-values that underflow to 0
    
    Interpretation:
        - 2 = p=0.01
        - 3 = p=0.001
        - 10 = p=1e-10
        - 300 = p < 1e-300 (underflow)
    """
    if p is None or np.isnan(p):
        return np.nan
    if p <= 0:
        return 300.0  # Cap at 10^-300
    if p >= 1:
        return 0.0
    return min(-np.log10(p), 300.0)


def format_pvalue(p):
    """
    Format p-value for display, handling very small values.
    
    Returns string like "p < 1e-300" or "p = 0.0234"
    """
    if p is None or np.isnan(p):
        return "N/A"
    if p <= 0:
        return "p < 1e-300"
    if p < 1e-300:
        return "p < 1e-300"
    if p < 0.0001:
        return f"p = {p:.2e}"
    return f"p = {p:.4f}"


def cramers_v(contingency_table):
    """
    Calculate Cramér's V effect size for chi-square test.
    
    Interpretation:
        - 0.1 = small effect
        - 0.3 = medium effect
        - 0.5 = large effect
    
    Args:
        contingency_table: numpy array (rows x cols)
    
    Returns:
        Cramér's V (0 to 1)
    """
    chi2 = stats.chi2_contingency(contingency_table)[0]
    n = contingency_table.sum()
    min_dim = min(contingency_table.shape) - 1
    if min_dim == 0 or n == 0:
        return 0.0
    return np.sqrt(chi2 / (n * min_dim))


def cramers_v_with_bias_correction(contingency_table):
    """
    Calculate bias-corrected Cramér's V (Bergsma, 2013).
    Better for smaller samples.
    """
    chi2 = stats.chi2_contingency(contingency_table)[0]
    n = contingency_table.sum()
    r, k = contingency_table.shape
    
    if n == 0 or min(r, k) <= 1:
        return 0.0
    
    phi2 = chi2 / n
    phi2_corrected = max(0, phi2 - ((r - 1) * (k - 1)) / (n - 1))
    r_corrected = r - ((r - 1) ** 2) / (n - 1)
    k_corrected = k - ((k - 1) ** 2) / (n - 1)
    
    denominator = min(r_corrected - 1, k_corrected - 1)
    if denominator <= 0:
        return 0.0
    
    return np.sqrt(phi2_corrected / denominator)


def epsilon_squared(h_stat, n):
    """
    Calculate epsilon squared effect size for Kruskal-Wallis test.
    
    Interpretation:
        - 0.01 = small effect
        - 0.06 = medium effect
        - 0.14 = large effect
    
    Args:
        h_stat: Kruskal-Wallis H statistic
        n: total sample size
    
    Returns:
        epsilon squared (0 to 1)
    """
    if n <= 1:
        return 0.0
    return h_stat / (n - 1)


def eta_squared_from_h(h_stat, n, k):
    """
    Calculate eta squared from Kruskal-Wallis H.
    
    Args:
        h_stat: Kruskal-Wallis H statistic
        n: total sample size
        k: number of groups
    
    Returns:
        eta squared (0 to 1)
    """
    if n <= k:
        return 0.0
    return (h_stat - k + 1) / (n - k)


def odds_ratio_2x2(a, b, c, d):
    """
    Calculate odds ratio for 2x2 contingency table.
    
    Table layout:
              Has_P2  No_P2
    Has_P1      a       b
    No_P1       c       d
    
    Returns:
        odds_ratio, 95% CI lower, 95% CI upper
    """
    # Add 0.5 to avoid division by zero (Haldane-Anscombe correction)
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    
    odds_ratio = (a * d) / (b * c)
    
    # Log odds ratio standard error
    se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
    
    # 95% CI
    log_or = np.log(odds_ratio)
    ci_lower = np.exp(log_or - 1.96 * se_log_or)
    ci_upper = np.exp(log_or + 1.96 * se_log_or)
    
    return odds_ratio, ci_lower, ci_upper


def cohens_h(p1, p2):
    """
    Calculate Cohen's h effect size for difference in proportions.
    
    Interpretation:
        - 0.2 = small effect
        - 0.5 = medium effect
        - 0.8 = large effect
    """
    phi1 = 2 * np.arcsin(np.sqrt(p1))
    phi2 = 2 * np.arcsin(np.sqrt(p2))
    return abs(phi1 - phi2)


def interpret_effect_size(value, metric='cramers_v'):
    """Return interpretation of effect size."""
    thresholds = {
        'cramers_v': [(0.1, 'negligible'), (0.3, 'small'), (0.5, 'medium'), (1.0, 'large')],
        'epsilon_sq': [(0.01, 'negligible'), (0.06, 'small'), (0.14, 'medium'), (1.0, 'large')],
        'cohens_h': [(0.2, 'negligible'), (0.5, 'small'), (0.8, 'medium'), (2.0, 'large')],
        'odds_ratio': [(1.5, 'negligible'), (2.5, 'small'), (4.0, 'medium'), (100.0, 'large')],
    }
    
    for threshold, label in thresholds.get(metric, [(1.0, 'unknown')]):
        if value < threshold:
            return label
    return 'large'


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data():
    """Load and prepare the annotated patient data."""
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} patients")
    print(f"Columns: {len(df.columns)}")
    
    # Verify cluster column exists
    if 'cluster' not in df.columns:
        raise ValueError("ERROR: 'cluster' column not found in input CSV!")
    
    # Get carriers only for most analyses (cluster != 3)
    carriers = df[df['n_pathways_affected'] > 0].copy()
    print(f"Carriers: {len(carriers)}")
    
    # Print cluster distribution
    print("\nCluster distribution (carriers only):")
    for cid in CARRIER_CLUSTERS:
        n = len(carriers[carriers['cluster'] == cid])
        label = CLUSTER_LABELS.get(cid, 'Unknown')
        print(f"  {cid} ({label}): {n}")
    
    # Print sample size warning
    if len(carriers) > 1000:
        print(f"\n⚠️  LARGE SAMPLE SIZE WARNING (n={len(carriers)})")
        print("   With this sample size, nearly all differences will be statistically")
        print("   significant. Focus on EFFECT SIZES rather than p-values.")
    
    return df, carriers


# =============================================================================
# 4.1: PATHWAY PREVALENCE BY SEVERITY
# =============================================================================
def analyze_pathway_prevalence(df, carriers):
    """Calculate pathway prevalence by Phase 1 cluster."""
    print("\n" + "="*70)
    print("4.1: PATHWAY PREVALENCE BY CLUSTER")
    print("="*70)
    
    results = []
    
    for cluster_id in CARRIER_CLUSTERS:
        cluster_carriers = carriers[carriers['cluster'] == cluster_id]
        n_carriers = len(cluster_carriers)
        
        if n_carriers == 0:
            continue
        
        cluster_label = CLUSTER_LABELS[cluster_id]
        row = {'cluster': cluster_id, 'cluster_label': cluster_label, 'n_carriers': n_carriers}
        
        for pathway in PATHWAYS:
            col = f'pathway_{pathway}'
            count = cluster_carriers[col].sum()
            pct = 100 * count / n_carriers if n_carriers > 0 else 0
            row[f'{pathway}_count'] = count
            row[f'{pathway}_pct'] = round(pct, 1)
        
        results.append(row)
        
        print(f"\n{cluster_label} (n={n_carriers}):")
        for pathway in PATHWAYS:
            print(f"  {pathway}: {row[f'{pathway}_count']} ({row[f'{pathway}_pct']:.1f}%)")
    
    # Create summary dataframe
    prevalence_df = pd.DataFrame(results)
    
    # Also create a pivot-friendly version
    pivot_data = []
    for row in results:
        for pathway in PATHWAYS:
            pivot_data.append({
                'cluster': row['cluster'],
                'cluster_label': row['cluster_label'],
                'pathway': pathway,
                'count': row[f'{pathway}_count'],
                'pct': row[f'{pathway}_pct'],
                'n_carriers': row['n_carriers']
            })
    
    pivot_df = pd.DataFrame(pivot_data)
    pivot_df.to_csv(f"{OUTPUT_DIR}/4.1_pathway_prevalence_by_cluster.csv", index=False)
    print(f"\nSaved: 4.1_pathway_prevalence_by_cluster.csv")
    
    return prevalence_df, pivot_df


# =============================================================================
# 4.2A: CO-OCCURRENCE MATRICES (BINARY) - WITH STATISTICS
# =============================================================================
def analyze_cooccurrence(carriers):
    """Calculate pathway co-occurrence matrices with statistical tests."""
    print("\n" + "="*70)
    print("4.2A: PATHWAY CO-OCCURRENCE MATRICES")
    print("="*70)
    
    all_results = {}
    
    # Overall co-occurrence
    print("\n--- Overall Co-occurrence ---")
    cooc_matrix = calculate_cooccurrence_matrix(carriers)
    all_results['Overall'] = cooc_matrix
    print(cooc_matrix.to_string())
    
    # Calculate statistical significance of co-occurrences
    print("\n--- Co-occurrence Statistics (Overall) ---")
    cooc_stats = calculate_cooccurrence_statistics(carriers)
    all_results['Overall_Stats'] = cooc_stats
    
    # Print top co-occurrences by effect size
    print("\nTop co-occurrences by odds ratio:")
    top_cooc = cooc_stats.nlargest(5, 'odds_ratio')
    for _, row in top_cooc.iterrows():
        print(f"  {row['pathway_1']} × {row['pathway_2']}: OR={row['odds_ratio']:.2f} "
              f"[{row['or_ci_lower']:.2f}-{row['or_ci_upper']:.2f}], "
              f"-log10(p)={row['neglog10_pvalue']:.1f}")
    
    # By cluster
    for cluster_id in CARRIER_CLUSTERS:
        cluster_carriers = carriers[carriers['cluster'] == cluster_id]
        if len(cluster_carriers) < 10:
            continue
        cluster_label = CLUSTER_LABELS[cluster_id]
        print(f"\n--- {cluster_label} (n={len(cluster_carriers)}) ---")
        cooc_matrix = calculate_cooccurrence_matrix(cluster_carriers)
        all_results[cluster_label] = cooc_matrix
        print(cooc_matrix.to_string())
    
    # Save all matrices to single file
    with pd.ExcelWriter(f"{OUTPUT_DIR}/4.2A_cooccurrence_matrices.xlsx") as writer:
        for name, matrix in all_results.items():
            if isinstance(matrix, pd.DataFrame):
                matrix.to_excel(writer, sheet_name=name[:31])  # Excel sheet name limit
    
    # Also save as CSV (overall only)
    all_results['Overall'].to_csv(f"{OUTPUT_DIR}/4.2A_cooccurrence_overall.csv")
    cooc_stats.to_csv(f"{OUTPUT_DIR}/4.2A_cooccurrence_statistics.csv", index=False)
    print(f"\nSaved: 4.2A_cooccurrence_matrices.xlsx, 4.2A_cooccurrence_overall.csv")
    print(f"Saved: 4.2A_cooccurrence_statistics.csv")
    
    return all_results


def calculate_cooccurrence_matrix(df):
    """Calculate pairwise co-occurrence counts."""
    n = len(PATHWAYS)
    matrix = np.zeros((n, n), dtype=int)
    
    for i, p1 in enumerate(PATHWAYS):
        col1 = f'pathway_{p1}'
        for j, p2 in enumerate(PATHWAYS):
            col2 = f'pathway_{p2}'
            # Count patients with both pathways
            count = ((df[col1] == 1) & (df[col2] == 1)).sum()
            matrix[i, j] = count
    
    return pd.DataFrame(matrix, index=PATHWAYS, columns=PATHWAYS)


def calculate_cooccurrence_statistics(df):
    """
    Calculate statistical significance of pathway co-occurrences.
    
    For each pair, calculates:
    - Fisher's exact test p-value
    - Odds ratio with 95% CI
    - Jaccard similarity
    """
    n_total = len(df)
    results = []
    
    for i, p1 in enumerate(PATHWAYS):
        for j, p2 in enumerate(PATHWAYS):
            if i >= j:  # Skip diagonal and lower triangle
                continue
            
            col1 = f'pathway_{p1}'
            col2 = f'pathway_{p2}'
            
            # Build 2x2 contingency table
            both = ((df[col1] == 1) & (df[col2] == 1)).sum()
            p1_only = ((df[col1] == 1) & (df[col2] == 0)).sum()
            p2_only = ((df[col1] == 0) & (df[col2] == 1)).sum()
            neither = ((df[col1] == 0) & (df[col2] == 0)).sum()
            
            # Fisher's exact test
            contingency = np.array([[both, p1_only], [p2_only, neither]])
            try:
                odds_ratio_fisher, fisher_p = stats.fisher_exact(contingency)
            except:
                odds_ratio_fisher, fisher_p = np.nan, np.nan
            
            # Calculate odds ratio with CI
            or_val, or_ci_low, or_ci_high = odds_ratio_2x2(both, p1_only, p2_only, neither)
            
            # Jaccard similarity
            union = both + p1_only + p2_only
            jaccard = both / union if union > 0 else 0
            
            # Expected co-occurrence under independence
            n1 = both + p1_only
            n2 = both + p2_only
            expected = (n1 * n2) / n_total if n_total > 0 else 0
            observed_expected_ratio = both / expected if expected > 0 else np.nan
            
            results.append({
                'pathway_1': p1,
                'pathway_2': p2,
                'cooccurrence_count': both,
                'expected_count': round(expected, 1),
                'obs_exp_ratio': round(observed_expected_ratio, 2) if not np.isnan(observed_expected_ratio) else None,
                'jaccard': round(jaccard, 3),
                'odds_ratio': round(or_val, 2),
                'or_ci_lower': round(or_ci_low, 2),
                'or_ci_upper': round(or_ci_high, 2),
                'fisher_pvalue': fisher_p,
                'neglog10_pvalue': round(safe_log10_pvalue(fisher_p), 1),
                'effect_interpretation': interpret_effect_size(or_val, 'odds_ratio'),
            })
    
    return pd.DataFrame(results)


# =============================================================================
# 4.2B: CORRELATION MATRICES (SCORES)
# =============================================================================
def analyze_correlations(carriers):
    """Calculate pathway score correlations (Spearman AND Pearson per Speaker 1 request)."""
    print("\n" + "="*70)
    print("4.2B: PATHWAY CORRELATION MATRICES")
    print("="*70)
    
    score_cols = [f'pathway_{p}_count' for p in PATHWAYS]
    
    all_results = {}
    
    # Overall Spearman correlations
    print("\n--- Overall Spearman Correlations ---")
    corr_spearman = carriers[score_cols].corr(method='spearman')
    corr_spearman.index = PATHWAYS
    corr_spearman.columns = PATHWAYS
    all_results['Overall_Spearman'] = corr_spearman
    print(corr_spearman.round(3).to_string())
    
    # Overall Pearson correlations (Speaker 1 requested for linear dose-dependent relationships)
    print("\n--- Overall Pearson Correlations ---")
    corr_pearson = carriers[score_cols].corr(method='pearson')
    corr_pearson.index = PATHWAYS
    corr_pearson.columns = PATHWAYS
    all_results['Overall_Pearson'] = corr_pearson
    print(corr_pearson.round(3).to_string())
    
    # Compare Spearman vs Pearson
    print("\n--- Spearman vs Pearson Comparison ---")
    diff = (corr_spearman - corr_pearson).abs()
    max_diff = diff.max().max()
    print(f"Max absolute difference: {max_diff:.3f}")
    if max_diff > 0.1:
        print("⚠️  Notable differences detected - relationships may be non-linear")
    else:
        print("✓ Correlations are similar - relationships appear linear")
    
    # Legacy reference for backward compatibility
    all_results['Overall'] = corr_spearman
    
    # By cluster (both Spearman and Pearson)
    for cluster_id in CARRIER_CLUSTERS:
        cluster_carriers = carriers[carriers['cluster'] == cluster_id]
        if len(cluster_carriers) < 30:
            continue
        cluster_label = CLUSTER_LABELS[cluster_id]
        print(f"\n--- {cluster_label} (n={len(cluster_carriers)}) ---")
        
        # Spearman
        corr_spearman = cluster_carriers[score_cols].corr(method='spearman')
        corr_spearman.index = PATHWAYS
        corr_spearman.columns = PATHWAYS
        all_results[f'{cluster_label}_Spearman'] = corr_spearman
        print("Spearman:")
        print(corr_spearman.round(3).to_string())
        
        # Pearson
        corr_pearson = cluster_carriers[score_cols].corr(method='pearson')
        corr_pearson.index = PATHWAYS
        corr_pearson.columns = PATHWAYS
        all_results[f'{cluster_label}_Pearson'] = corr_pearson
    
    # Save all matrices
    with pd.ExcelWriter(f"{OUTPUT_DIR}/4.2B_correlation_matrices.xlsx") as writer:
        for name, matrix in all_results.items():
            # Excel sheet names limited to 31 chars
            sheet_name = name[:31]
            matrix.to_excel(writer, sheet_name=sheet_name)
    
    all_results['Overall_Spearman'].to_csv(f"{OUTPUT_DIR}/4.2B_correlation_spearman.csv")
    all_results['Overall_Pearson'].to_csv(f"{OUTPUT_DIR}/4.2B_correlation_pearson.csv")
    print(f"\nSaved: 4.2B_correlation_matrices.xlsx")
    print(f"Saved: 4.2B_correlation_spearman.csv, 4.2B_correlation_pearson.csv")
    
    return all_results


# =============================================================================
# 4.3: FREQUENCY VS INTENSITY COMPARISON
# =============================================================================
def analyze_frequency_vs_intensity(carriers):
    """Compare pathway frequency (binary) vs intensity (score) by cluster."""
    print("\n" + "="*70)
    print("4.3: FREQUENCY VS INTENSITY ANALYSIS")
    print("="*70)
    
    results = []
    
    for cluster_id in CARRIER_CLUSTERS:
        cluster_carriers = carriers[carriers['cluster'] == cluster_id]
        n = len(cluster_carriers)
        if n == 0:
            continue
        
        cluster_label = CLUSTER_LABELS[cluster_id]
        print(f"\n--- {cluster_label} (n={n}) ---")
        
        for pathway in PATHWAYS:
            binary_col = f'pathway_{pathway}'
            score_col = f'pathway_{pathway}_count'
            
            # Frequency: % of carriers with this pathway
            freq = 100 * cluster_carriers[binary_col].sum() / n
            
            # Intensity: mean score among those WITH the pathway
            affected = cluster_carriers[cluster_carriers[binary_col] == 1]
            intensity = affected[score_col].mean() if len(affected) > 0 else 0
            
            # Also calculate overall mean score (including zeros)
            overall_intensity = cluster_carriers[score_col].mean()
            
            results.append({
                'cluster': cluster_id,
                'cluster_label': cluster_label,
                'pathway': pathway,
                'n_carriers': n,
                'frequency_pct': round(freq, 1),
                'mean_count_affected': round(intensity, 2),
                'mean_count_overall': round(overall_intensity, 2),
            })
            
            print(f"  {pathway}: {freq:.1f}% freq, {intensity:.2f} mean score (affected)")
    
    freq_int_df = pd.DataFrame(results)
    freq_int_df.to_csv(f"{OUTPUT_DIR}/4.3_frequency_vs_intensity.csv", index=False)
    print(f"\nSaved: 4.3_frequency_vs_intensity.csv")
    
    return freq_int_df


# =============================================================================
# 4.4: NETWORK VISUALIZATION DATA
# =============================================================================
def generate_network_data(carriers, cooc_matrices):
    """Generate network visualization data for pathway relationships."""
    print("\n" + "="*70)
    print("4.4: NETWORK VISUALIZATION DATA")
    print("="*70)
    
    n_total = len(carriers)
    
    # Create nodes (pathways)
    nodes = []
    for pathway in PATHWAYS:
        col = f'pathway_{pathway}'
        count = carriers[col].sum()
        nodes.append({
            'id': pathway,
            'label': pathway.replace('_', ' '),
            'size': int(count),
            'pct': round(100 * count / n_total, 1)
        })
    
    # Create edges (co-occurrences) with effect sizes
    edges = []
    cooc_overall = cooc_matrices['Overall']
    
    for i, p1 in enumerate(PATHWAYS):
        for j, p2 in enumerate(PATHWAYS):
            if i < j:  # Upper triangle only
                count = cooc_overall.loc[p1, p2]
                if count > 0:
                    # Calculate Jaccard similarity
                    n1 = cooc_overall.loc[p1, p1]
                    n2 = cooc_overall.loc[p2, p2]
                    jaccard = count / (n1 + n2 - count) if (n1 + n2 - count) > 0 else 0
                    
                    # Calculate odds ratio
                    col1 = f'pathway_{p1}'
                    col2 = f'pathway_{p2}'
                    both = count
                    p1_only = n1 - both
                    p2_only = n2 - both
                    neither = n_total - n1 - n2 + both
                    
                    or_val, or_ci_low, or_ci_high = odds_ratio_2x2(both, p1_only, p2_only, neither)
                    
                    # Expected under independence
                    expected = (n1 * n2) / n_total if n_total > 0 else 0
                    
                    edges.append({
                        'source': p1,
                        'target': p2,
                        'weight': int(count),
                        'expected': round(expected, 1),
                        'jaccard': round(jaccard, 3),
                        'odds_ratio': round(or_val, 2),
                        'or_ci_lower': round(or_ci_low, 2),
                        'or_ci_upper': round(or_ci_high, 2),
                        'effect_size': interpret_effect_size(or_val, 'odds_ratio'),
                    })
    
    network_data = {'nodes': nodes, 'edges': edges}
    
    # Get Pearson correlations for quadrant assignment
    score_cols = [f'pathway_{p}_count' for p in PATHWAYS]
    corr_pearson = carriers[score_cols].corr(method='pearson')
    corr_pearson.index = PATHWAYS
    corr_pearson.columns = PATHWAYS
    
    # Add quadrant metadata for visualization (matches Speaker 1's scatter plot)
    # X-axis: Co-occurrence % (of smaller pathway)
    # Y-axis: Pearson correlation
    network_data['quadrant_config'] = {
        'x_axis': {
            'label': 'Co-occurrence (%)',
            'description': 'Percentage of smaller pathway that co-occurs with larger'
        },
        'y_axis': {
            'label': 'Pearson Correlation',
            'description': 'Linear correlation between pathway counts'
        },
        'categories': {
            'dose_dependent': {
                'label': 'Dose-dependent (cascading failure)',
                'color': 'green',
                'criteria': 'co-occ >50%, r >0.5',
                'interpretation': 'Strong linear relationship - pathway burden compounds'
            },
            'threshold_effect': {
                'label': 'Threshold effect (binary disruption)',
                'color': 'orange',
                'criteria': 'co-occ >50%, |r| <0.3',
                'interpretation': 'High co-occurrence but weak correlation - presence matters, not dose'
            },
            'distinct_subtypes': {
                'label': 'Distinct subtypes',
                'color': 'red',
                'criteria': 'co-occ <30%, r <-0.3',
                'interpretation': 'Negative correlation - mutually exclusive pathway patterns'
            },
            'other': {
                'label': 'Other patterns',
                'color': 'gray',
                'criteria': 'Does not meet above criteria'
            }
        }
    }
    
    # Assign quadrant category to each edge
    for edge in edges:
        src = edge['source']
        tgt = edge['target']
        
        # Get co-occurrence percentage (relative to smaller pathway)
        src_n = next((n['size'] for n in nodes if n['id'] == src), 0)
        tgt_n = next((n['size'] for n in nodes if n['id'] == tgt), 0)
        min_n = min(src_n, tgt_n)
        cooc_pct = 100 * edge['weight'] / min_n if min_n > 0 else 0
        edge['cooc_pct'] = round(cooc_pct, 1)
        
        # Get Pearson correlation
        r = corr_pearson.loc[src, tgt]
        edge['pearson_r'] = round(r, 3)
        
        # Assign category based on criteria from the plot
        if cooc_pct > 50 and r > 0.5:
            edge['category'] = 'dose_dependent'
            edge['color'] = 'green'
        elif cooc_pct > 50 and abs(r) < 0.3:
            edge['category'] = 'threshold_effect'
            edge['color'] = 'orange'
        elif cooc_pct < 30 and r < -0.3:
            edge['category'] = 'distinct_subtypes'
            edge['color'] = 'red'
        else:
            edge['category'] = 'other'
            edge['color'] = 'gray'
    
    # Save as JSON for visualization tools
    def convert_types(obj):
        if isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(i) for i in obj]
        elif isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64)):
            return float(obj)
        return obj
    
    with open(f"{OUTPUT_DIR}/4.4_network_data.json", 'w') as f:
        json.dump(convert_types(network_data), f, indent=2)
    
    # Also save as CSV
    pd.DataFrame(nodes).to_csv(f"{OUTPUT_DIR}/4.4_network_nodes.csv", index=False)
    pd.DataFrame(edges).to_csv(f"{OUTPUT_DIR}/4.4_network_edges.csv", index=False)
    
    print(f"\nNodes: {len(nodes)}")
    print(f"Edges: {len(edges)}")
    print(f"\nTop 5 co-occurrences by EFFECT SIZE (odds ratio):")
    sorted_edges = sorted(edges, key=lambda x: x['odds_ratio'], reverse=True)[:5]
    for e in sorted_edges:
        print(f"  {e['source']} <-> {e['target']}: n={e['weight']}, OR={e['odds_ratio']:.2f} ({e['effect_size']})")
    
    # Print category distribution (matches scatter plot quadrants)
    print(f"\n--- Category Distribution (Co-occurrence % vs Pearson r) ---")
    categories = {'dose_dependent': [], 'threshold_effect': [], 'distinct_subtypes': [], 'other': []}
    for e in edges:
        categories[e['category']].append(e)
    
    print(f"\n🟢 Dose-dependent (co-occ >50%, r >0.5) - Cascading failure:")
    for e in categories['dose_dependent']:
        print(f"     {e['source']} × {e['target']}: {e['cooc_pct']:.0f}% co-occ, r={e['pearson_r']:.2f}")
    if not categories['dose_dependent']:
        print(f"     (none)")
    
    print(f"\n🟠 Threshold effect (co-occ >50%, |r| <0.3) - Binary disruption:")
    for e in categories['threshold_effect']:
        print(f"     {e['source']} × {e['target']}: {e['cooc_pct']:.0f}% co-occ, r={e['pearson_r']:.2f}")
    if not categories['threshold_effect']:
        print(f"     (none)")
    
    print(f"\n🔴 Distinct subtypes (co-occ <30%, r <-0.3) - Independent:")
    for e in categories['distinct_subtypes']:
        print(f"     {e['source']} × {e['target']}: {e['cooc_pct']:.0f}% co-occ, r={e['pearson_r']:.2f}")
    if not categories['distinct_subtypes']:
        print(f"     (none)")
    
    print(f"\n⚪ Other patterns: {len(categories['other'])} pairs")
    
    print(f"\nSaved: 4.4_network_data.json, 4.4_network_nodes.csv, 4.4_network_edges.csv")
    
    return network_data


# =============================================================================
# 4.5: CROSS-CLUSTER COMPARISON (WITH EFFECT SIZES)
# =============================================================================
def analyze_cross_cluster(carriers):
    """Statistical comparison of pathways across Phase 1 clusters with effect sizes."""
    print("\n" + "="*70)
    print("4.5: CROSS-CLUSTER STATISTICAL COMPARISON")
    print("(With effect sizes - these matter more than p-values!)")
    print("="*70)
    
    results = []
    n_total = len(carriers)
    
    for pathway in PATHWAYS:
        binary_col = f'pathway_{pathway}'
        score_col = f'pathway_{pathway}_count'
        
        # Get data by cluster
        cluster_data = {}
        for cluster_id in CARRIER_CLUSTERS:
            cluster_df = carriers[carriers['cluster'] == cluster_id]
            cluster_label = CLUSTER_LABELS[cluster_id]
            cluster_data[cluster_id] = {
                'label': cluster_label,
                'n': len(cluster_df),
                'has_pathway': cluster_df[binary_col].sum(),
                'no_pathway': len(cluster_df) - cluster_df[binary_col].sum(),
                'scores': cluster_df[score_col].values,
                'prevalence': cluster_df[binary_col].mean() if len(cluster_df) > 0 else 0,
            }
        
        # Chi-square test for binary (across all clusters)
        contingency_data = []
        for cluster_id in CARRIER_CLUSTERS:
            contingency_data.append([
                cluster_data[cluster_id]['has_pathway'],
                cluster_data[cluster_id]['no_pathway']
            ])
        
        contingency = np.array(contingency_data)
        
        if contingency.min() >= 5:  # Chi-square assumption
            chi2, chi_p, dof, expected = stats.chi2_contingency(contingency)
            cramers = cramers_v(contingency)
            cramers_corrected = cramers_v_with_bias_correction(contingency)
        else:
            chi2, chi_p = np.nan, np.nan
            cramers = np.nan
            cramers_corrected = np.nan
        
        # Kruskal-Wallis test for scores
        groups = []
        for cluster_id in CARRIER_CLUSTERS:
            scores = cluster_data[cluster_id]['scores']
            if len(scores) > 0:
                groups.append(scores)
        
        if len(groups) >= 2 and all(len(g) > 0 for g in groups):
            h_stat, kw_p = stats.kruskal(*groups)
            epsilon_sq = epsilon_squared(h_stat, n_total)
            eta_sq = eta_squared_from_h(h_stat, n_total, len(groups))
        else:
            h_stat, kw_p = np.nan, np.nan
            epsilon_sq = np.nan
            eta_sq = np.nan
        
        # Pairwise effect sizes (Cluster 0 vs Cluster 2)
        p_cluster0 = cluster_data[0]['prevalence']
        p_cluster2 = cluster_data[2]['prevalence']
        cohens_h_0_vs_2 = cohens_h(p_cluster0, p_cluster2) if p_cluster0 > 0 or p_cluster2 > 0 else 0
        
        # Prevalence values per cluster
        prev_c0 = 100 * cluster_data[0]['prevalence']
        prev_c1 = 100 * cluster_data[1]['prevalence']
        prev_c2 = 100 * cluster_data[2]['prevalence']
        prev_c4 = 100 * cluster_data[4]['prevalence']
        
        # Trend direction (comparing Cluster 0 mild to Cluster 2 severe)
        if prev_c2 > prev_c0 + 5:
            trend = "↑ increases mild→severe"
        elif prev_c2 < prev_c0 - 5:
            trend = "↓ decreases mild→severe"
        else:
            trend = "→ no clear trend"
        
        results.append({
            'pathway': pathway,
            'n_cluster_0': cluster_data[0]['n'],
            'n_cluster_1': cluster_data[1]['n'],
            'n_cluster_2': cluster_data[2]['n'],
            'n_cluster_4': cluster_data[4]['n'],
            'prevalence_cluster_0_pct': round(prev_c0, 1),
            'prevalence_cluster_1_pct': round(prev_c1, 1),
            'prevalence_cluster_2_pct': round(prev_c2, 1),
            'prevalence_cluster_4_pct': round(prev_c4, 1),
            'trend': trend,
            # Chi-square results
            'chi2_statistic': round(chi2, 2) if not np.isnan(chi2) else None,
            'chi2_neglog10p': round(safe_log10_pvalue(chi_p), 1) if not np.isnan(chi_p) else None,
            'cramers_v': round(cramers, 3) if not np.isnan(cramers) else None,
            'cramers_v_corrected': round(cramers_corrected, 3) if not np.isnan(cramers_corrected) else None,
            'chi2_effect_interpretation': interpret_effect_size(cramers, 'cramers_v') if not np.isnan(cramers) else None,
            # Kruskal-Wallis results  
            'kruskal_H': round(h_stat, 2) if not np.isnan(h_stat) else None,
            'kruskal_neglog10p': round(safe_log10_pvalue(kw_p), 1) if not np.isnan(kw_p) else None,
            'epsilon_squared': round(epsilon_sq, 4) if not np.isnan(epsilon_sq) else None,
            'eta_squared': round(eta_sq, 4) if not np.isnan(eta_sq) else None,
            'kruskal_effect_interpretation': interpret_effect_size(epsilon_sq, 'epsilon_sq') if not np.isnan(epsilon_sq) else None,
            # Pairwise effect size
            'cohens_h_c0_vs_c2': round(cohens_h_0_vs_2, 3),
            'cohens_h_interpretation': interpret_effect_size(cohens_h_0_vs_2, 'cohens_h'),
        })
        
        # Print summary
        effect_chi = interpret_effect_size(cramers, 'cramers_v') if not np.isnan(cramers) else "N/A"
        effect_kw = interpret_effect_size(epsilon_sq, 'epsilon_sq') if not np.isnan(epsilon_sq) else "N/A"
        
        print(f"\n{pathway}:")
        print(f"  Prevalence: C0={prev_c0:.1f}%, C1={prev_c1:.1f}%, C2={prev_c2:.1f}%, C4={prev_c4:.1f}% {trend}")
        if not np.isnan(chi2):
            print(f"  Chi-square: χ²={chi2:.2f}, -log10(p)={safe_log10_pvalue(chi_p):.1f}, Cramér's V={cramers:.3f} ({effect_chi})")
        else:
            print(f"  Chi-square: N/A")
        if not np.isnan(h_stat):
            print(f"  Kruskal-Wallis: H={h_stat:.2f}, -log10(p)={safe_log10_pvalue(kw_p):.1f}, ε²={epsilon_sq:.4f} ({effect_kw})")
        else:
            print(f"  Kruskal-Wallis: N/A")
        print(f"  Cohen's h (C0 vs C2): {cohens_h_0_vs_2:.3f} ({interpret_effect_size(cohens_h_0_vs_2, 'cohens_h')})")
    
    comparison_df = pd.DataFrame(results)
    comparison_df.to_csv(f"{OUTPUT_DIR}/4.5_cross_cluster_comparison.csv", index=False)
    print(f"\nSaved: 4.5_cross_cluster_comparison.csv")
    
    # Print summary of clinically meaningful effects
    print("\n" + "-"*50)
    print("SUMMARY: Pathways with medium/large effect sizes:")
    print("-"*50)
    
    meaningful = comparison_df[
        (comparison_df['chi2_effect_interpretation'].isin(['medium', 'large'])) |
        (comparison_df['kruskal_effect_interpretation'].isin(['medium', 'large'])) |
        (comparison_df['cohens_h_interpretation'].isin(['medium', 'large']))
    ]
    
    if len(meaningful) > 0:
        for _, row in meaningful.iterrows():
            print(f"  {row['pathway']}: {row['trend']}")
            print(f"    Cramér's V = {row['cramers_v']} ({row['chi2_effect_interpretation']})")
    else:
        print("  No pathways showed medium or large effect sizes across severity.")
        print("  All differences, while statistically significant, are clinically small.")
    
    return comparison_df


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("="*70)
    print("STAGE 4: PATHWAY CO-OCCURRENCE AND CORRELATION ANALYSIS")
    print("(CLUSTER VERSION - using Phase 1 clusters with effect sizes)")
    print("="*70)
    
    # Load data
    df, carriers = load_data()
    
    # Run all analyses
    prevalence_df, pivot_df = analyze_pathway_prevalence(df, carriers)
    cooc_matrices = analyze_cooccurrence(carriers)
    corr_matrices = analyze_correlations(carriers)
    freq_int_df = analyze_frequency_vs_intensity(carriers)
    network_data = generate_network_data(carriers, cooc_matrices)
    comparison_df = analyze_cross_cluster(carriers)
    
    # Summary
    print("\n" + "="*70)
    print("STAGE 4 COMPLETE - OUTPUT FILES")
    print("="*70)
    print(f"""
    4.1_pathway_prevalence_by_cluster.csv
    4.2A_cooccurrence_matrices.xlsx
    4.2A_cooccurrence_overall.csv
    4.2A_cooccurrence_statistics.csv  ← includes odds ratios & effect sizes
    4.2B_correlation_matrices.xlsx    ← includes BOTH Spearman and Pearson
    4.2B_correlation_spearman.csv
    4.2B_correlation_pearson.csv      ← NEW: per Speaker 1 request
    4.3_frequency_vs_intensity.csv
    4.4_network_data.json             ← UPDATED: includes quadrant labels
    4.4_network_nodes.csv
    4.4_network_edges.csv
    4.5_cross_cluster_comparison.csv  ← includes effect sizes
    """)
    
    print("="*70)
    print("KEY CHANGES IN THIS PATCHED VERSION:")
    print("="*70)
    print("""
    1. P-values reported as -log10(p) to handle underflow
       - Value of 10 = p=1e-10
       - Value of 300 = p < 1e-300 (underflow)
    
    2. Effect sizes added for clinical interpretation:
       - Cramér's V: effect size for chi-square (0.1=small, 0.3=medium, 0.5=large)
       - Epsilon squared: effect size for Kruskal-Wallis (0.01=small, 0.06=medium, 0.14=large)
       - Odds ratio: effect size for co-occurrence (>2.5=small, >4=medium)
       - Cohen's h: effect size for prevalence differences (0.2=small, 0.5=medium, 0.8=large)
    
    3. Network data now includes odds ratios for edge weights
    
    4. Effect interpretations provided (negligible/small/medium/large)
    
    REMEMBER: With large sample sizes, focus on EFFECT SIZES, not p-values!
    """)
    
    return {
        'prevalence': prevalence_df,
        'cooccurrence': cooc_matrices,
        'correlations': corr_matrices,
        'freq_intensity': freq_int_df,
        'network': network_data,
        'cross_severity': comparison_df
    }


if __name__ == "__main__":
    main()