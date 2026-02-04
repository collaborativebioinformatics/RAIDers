# Stage 4: Pathway Co-occurrence and Correlation Analysis
## Methodology Documentation (Cluster-Based Version)

**Analysis Date:** February 2025  
**Input:** `patients_with_pathways_weighted_clusters.csv` (N=15,000 patients, 6,043 carriers)  
**Clustering:** Phase 1 Federated K-means (4 carrier clusters)

---

## 1. Cluster Definitions

| Cluster | Label | N | Severity | Dominant Pathways |
|---------|-------|---|----------|-------------------|
| 0 | Cluster_0_Mild | 245 | Mild | Vesicle (45.7%), DNA Damage (42.0%) |
| 1 | Cluster_1_Moderate_SOD1 | 3,412 | Moderate | Proteostasis (64.6%), Excitotoxicity (50.4%) |
| 2 | Cluster_2_Severe | 1,767 | Severe | Proteostasis (72.9%), RNA Metabolism (26.9%) |
| 4 | Cluster_4_Moderate_ALS2 | 619 | Moderate | Vesicle (32.3%), Proteostasis (24.7%) |

**Key Insight:** Clusters 1 and 4 are both "Moderate" severity but represent genetically distinct subtypes:
- **C1 (SOD1-dominant):** High Proteostasis (64.6%) + Excitotoxicity (50.4%)
- **C4 (ALS2-dominant):** Higher Cytoskeletal (12.0% vs 4.3%), lower Proteostasis (24.7%)

---

## 2. Statistical Methods

### 2.1 Correlation Analysis (Section 4.2B)

**Both Spearman and Pearson correlations computed per Speaker 1 request:**

| Pathway Pair | Spearman ρ | Pearson r | Δ | Interpretation |
|--------------|------------|-----------|---|----------------|
| Mitochondrial × Excitotoxicity | 0.704 | 0.687 | 0.017 | Strong positive (linear) |
| Proteostasis × Excitotoxicity | 0.437 | 0.430 | 0.007 | Moderate positive |
| Vesicle × DNA_Damage | 0.364 | 0.375 | 0.011 | Moderate positive |
| Proteostasis × Vesicle | -0.442 | -0.408 | 0.034 | Moderate negative |
| Proteostasis × RNA_Metabolism | -0.258 | -0.240 | 0.018 | Weak negative |

**Max |Spearman - Pearson| = 0.034** → Relationships are approximately linear; both methods yield similar conclusions.

**Rationale for both methods:**
- **Spearman:** Robust to outliers, captures monotonic relationships
- **Pearson:** Detects linear dose-dependent relationships (Speaker 1's biological hypothesis)

### 2.2 Co-occurrence Analysis (Section 4.2A)

**Top associations by Odds Ratio:**

| Pathway Pair | Co-occurrence | OR | 95% CI | Effect |
|--------------|---------------|-----|--------|--------|
| Mitochondrial × Excitotoxicity | 1,769 | 38.78 | [33.23, 45.27] | Large |
| Proteostasis × Excitotoxicity | 1,821 | 9.63 | [8.28, 11.21] | Large |
| Vesicle × DNA_Damage | 627 | 6.94 | [5.98, 8.06] | Large |
| Proteostasis × Mitochondrial | 1,844 | 4.00 | [3.55, 4.50] | Medium |
| RNA_Metabolism × Mitochondrial | 479 | 2.79 | [2.40, 3.25] | Medium |

### 2.3 Cross-Cluster Comparison (Section 4.5)

**Effect sizes for pathway discrimination across clusters:**

| Pathway | ε² (Kruskal-Wallis) | Interpretation | Trend |
|---------|---------------------|----------------|-------|
| Excitotoxicity | 0.1614 | **Large** | ↑ increases with severity |
| Proteostasis | 0.1402 | **Large** | ↑ increases with severity |
| RNA_Metabolism | 0.0724 | Medium | ↑ increases with severity |
| Mitochondrial | 0.0580 | Small | ↑ increases with severity |
| Vesicle_Trafficking | 0.0473 | Small | ↓ decreases with severity |
| DNA_Damage | 0.0310 | Small | ↓ decreases with severity |
| Cytoskeletal | 0.0165 | Small | → no clear trend |

**Effect size interpretation (epsilon squared):**
- < 0.01: Negligible
- 0.01-0.06: Small
- 0.06-0.14: Medium
- > 0.14: Large

---

## 3. Network Visualization (Section 4.4)

### 3.1 Quadrant Classification

**Axes (matching scatter plot):**
- X-axis: Co-occurrence % (relative to smaller pathway)
- Y-axis: Pearson Correlation

| Category | Color | Criteria | Interpretation |
|----------|-------|----------|----------------|
| **Dose-dependent** | 🟢 Green | co-occ >50%, r >0.5 | Cascading failure - pathway burden compounds |
| **Threshold effect** | 🟠 Orange | co-occ >50%, \|r\| <0.3 | Binary disruption - presence matters, not dose |
| **Distinct subtypes** | 🔴 Red | co-occ <30%, r <-0.3 | Independent/mutually exclusive patterns |
| **Other** | ⚪ Gray | Does not meet above | Mixed or transitional patterns |

### 3.2 Category Assignments

**🟢 Dose-dependent (Cascading Failure):**
| Pair | Co-occ % | Pearson r |
|------|----------|-----------|
| Mitochondrial × Excitotoxicity | 86% | 0.69 |

*Interpretation:* Strong linear dose-response. As mitochondrial dysfunction increases, excitotoxicity increases proportionally. Suggests shared upstream drivers (e.g., SOD1 affecting both pathways).

**🟠 Threshold Effect (Binary Disruption):**
| Pair | Co-occ % | Pearson r |
|------|----------|-----------|
| RNA_Metabolism × Mitochondrial | 60% | 0.21 |

*Interpretation:* High co-occurrence but weak correlation. Both pathways are frequently affected together, but the *intensity* of one doesn't predict the other. Presence/absence matters more than severity.

**🔴 Distinct Subtypes:**
| Pair | Co-occ % | Pearson r |
|------|----------|-----------|
| Proteostasis × Vesicle_Trafficking | 24% | -0.41 |

*Interpretation:* Negative correlation with low overlap. Patients tend toward EITHER proteostasis-dominant OR vesicle-dominant phenotypes. Supports distinct molecular subtypes within ALS.

---

## 4. Key Findings

### 4.1 Cluster Validation (Addresses Peer Review Feedback)

Phase 1 clusters show distinct pathway signatures, validating biological substructure:

1. **C0 (Mild):** Isolated to Vesicle/DNA_Damage pathways (no Proteostasis/Excitotoxicity)
2. **C1 (SOD1-Moderate):** Core pathways dominant (Proteostasis + Excitotoxicity + Mitochondrial)
3. **C2 (Severe):** Highest Proteostasis (72.9%) + elevated RNA_Metabolism (26.9%)
4. **C4 (ALS2-Moderate):** Mixed pattern, highest Cytoskeletal involvement (12.0%)

### 4.2 Mitochondrial-Excitotoxicity Axis

Strongest association in dataset (OR=38.78, ρ=0.704):
- **Biological mechanism:** SOD1 mutations cause mitochondrial ROS accumulation → enhanced caspase-3 cleavage of EAAT2 → impaired glutamate clearance → excitotoxicity
- **Clinical implication:** Combined targeting may be more effective than single-pathway approaches

### 4.3 Negative Correlations Suggest Distinct Subtypes

Proteostasis shows negative correlation with Vesicle_Trafficking (ρ=-0.442):
- Patients tend toward EITHER proteostasis-dominant OR vesicle-dominant phenotypes
- Supports existence of distinct molecular subtypes within ALS

---

## 5. Output Files

| File | Description |
|------|-------------|
| `4.1_pathway_prevalence_by_cluster.csv` | Pathway prevalence per cluster |
| `4.2A_cooccurrence_matrices.xlsx` | Co-occurrence matrices (overall + per cluster) |
| `4.2A_cooccurrence_statistics.csv` | Odds ratios, Fisher's exact, effect sizes |
| `4.2B_correlation_matrices.xlsx` | Spearman AND Pearson correlations |
| `4.2B_correlation_spearman.csv` | Overall Spearman matrix |
| `4.2B_correlation_pearson.csv` | Overall Pearson matrix |
| `4.3_frequency_vs_intensity.csv` | Frequency vs mean gene count |
| `4.4_network_data.json` | Network viz with quadrant labels |
| `4.5_cross_cluster_comparison.csv` | Cross-cluster statistics with effect sizes |

---

## 6. Methodological Notes

### Sample Size Consideration
With N=6,043 carriers, nearly all tests achieve statistical significance. **Effect sizes are more informative than p-values** for this analysis.

### Correlation Method Choice
Both Spearman and Pearson included because:
- Spearman: Standard for non-normal count data
- Pearson: Requested by Speaker 1 for linear dose-dependent interpretation
- Results are concordant (max difference = 0.034)

### Cluster Assignment
Clusters derived from Phase 1 Federated K-means on severity scores. Assignment is deterministic and pre-computed in input CSV.

---

*Generated by stage4_pathway_analysis_clusters.py*
