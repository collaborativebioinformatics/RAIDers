# Methods

## 2.1 Study Design and Overview

This study employed a three-phase analytical pipeline to investigate molecular pathway heterogeneity in amyotrophic lateral sclerosis (ALS) using synthetic patient populations derived from clinically validated pathogenic variants (Figure 1). The methodology comprised: (1) synthetic patient generation with federated clustering to identify disease severity strata, (2) mechanistic pathway annotation based on curated gene-pathway relationships, and (3) comprehensive pathway co-occurrence and correlation analysis stratified by cluster assignment.

---

## 2.2 Data Sources

### 2.2.1 ClinVar Variant Database

Pathogenic ALS variants were extracted from the National Center for Biotechnology Information (NCBI) ClinVar database. The curated dataset (`clinvar.cleaned.csv`) was filtered using the following inclusion criteria:

**Inclusion Criteria:**
- Clinical significance classified as "Pathogenic" or "Likely Pathogenic" per American College of Medical Genetics and Genomics (ACMG) guidelines
- Disease association with amyotrophic lateral sclerosis or related motor neuron disease phenotypes
- Minimum review status of "criteria provided, single submitter"

**Resulting Dataset Characteristics:**

| Parameter | Value |
|-----------|-------|
| Total pathogenic variants | ~450 |
| Unique genes represented | 34 |
| Variant types included | SNV, indel, repeat expansion |

### 2.2.2 Variant Record Structure

Each ClinVar record provided the following biological ground truth for simulation:

| Field | Description | Example |
|-------|-------------|---------|
| Gene Association | HGNC-approved gene symbol | *SOD1*, *TARDBP*, *C9orf72* |
| Clinical Significance | ACMG/AMP pathogenicity classification | Pathogenic, Likely pathogenic |
| Molecular Consequence | Predicted functional impact | Missense, nonsense, frameshift |
| Genomic Coordinates | GRCh38 chromosome position | chr21:31659666 |

---

## 2.3 Synthetic Patient Population Generation

### 2.3.1 Population Architecture

Synthetic patient cohorts were generated to represent five continental superpopulations, following the 1000 Genomes Project nomenclature:

| Code | Population | Simulated n |
|------|------------|-------------|
| AFR | African/African American | 3,000 |
| AMR | Admixed American | 3,000 |
| EAS | East Asian | 3,000 |
| EUR | European | 3,000 |
| SAS | South Asian | 3,000 |
| **Total** | | **15,000** |

### 2.3.2 Variant Assignment Model

Pathogenic variants were assigned to synthetic patients using a stochastic model incorporating population-specific allele frequencies and gene-level penetrance estimates. For each patient *i* and variant *v*:

$$
P(\text{carries } v_j) = f_{v,pop} \times \pi_{gene(v)}
$$

Where:
- $f_{v,pop}$ = population-specific allele frequency from gnomAD v3.1
- $\pi_{gene(v)}$ = penetrance coefficient for the variant's associated gene

### 2.3.3 Severity Score Computation

Each variant carrier was assigned a composite severity score reflecting the cumulative pathogenic burden. The severity score for patient *i* carrying variants $\{v_1, v_2, \ldots, v_n\}$ was computed as:

$$
S_i = \sum_{j=1}^{n_i} w_{gene(v_j)} \times w_{consequence(v_j)} \times w_{pathogenicity(v_j)}
$$

**Component Weights:**

**Gene Weight ($w_{gene}$):** Literature-derived penetrance estimates

| Gene | Weight | Justification |
|------|--------|---------------|
| *SOD1* | 1.0 | High penetrance, well-characterized |
| *C9orf72* | 1.0 | Most common familial ALS cause |
| *TARDBP* | 0.9 | High penetrance |
| *FUS* | 0.9 | Aggressive juvenile-onset forms |
| Other genes | 0.5–0.8 | Variable penetrance |

**Molecular Consequence Weight ($w_{consequence}$):**

$$
w_{consequence} = \begin{cases}
1.0 & \text{if loss-of-function (nonsense, frameshift)} \\
0.9 & \text{if canonical splice site} \\
0.7 & \text{if missense with CADD} \geq 25 \\
0.5 & \text{if missense with } 15 \leq \text{CADD} < 25
\end{cases}
$$

**Pathogenicity Confidence Weight ($w_{pathogenicity}$):**

$$
w_{pathogenicity} = \begin{cases}
1.0 & \text{if ClinVar = "Pathogenic"} \\
0.8 & \text{if ClinVar = "Likely pathogenic"}
\end{cases}
$$

### 2.3.4 Severity Categorization

Continuous severity scores were discretized into clinical severity categories:

| Category | Severity Score Range | Clinical Interpretation |
|----------|---------------------|------------------------|
| Healthy | $S = 0$ | No pathogenic variants detected |
| Mild | $0 < S \leq 5.5$ | Single low-penetrance variant |
| Moderate | $5.5 < S \leq 8.0$ | Multiple variants or single high-impact |
| Severe | $S > 8.0$ | Multiple high-impact variants |

### 2.3.5 Predicted Progression Assignment

Disease progression rate was modeled as a probabilistic function of severity score:

$$
P(\text{progression} = k \mid S) = \text{softmax}\left(\beta_k \cdot S + \alpha_k\right)
$$

Where $k \in \{\text{Slow}, \text{Moderate}, \text{Fast}\}$ and parameters were calibrated to clinical literature on ALS progression rates.

---

## 2.4 Phase 1: Federated K-Means Clustering

### 2.4.1 Rationale

Federated learning was employed to simulate a privacy-preserving multi-institutional collaboration where patient-level data remains decentralized at each population node. This architecture enables identification of disease subgroups without sharing raw genetic data.

### 2.4.2 Feature Engineering

For each carrier patient, a feature vector $\mathbf{x}_i \in \mathbb{R}^d$ was constructed:

$$
\mathbf{x}_i = \left[ n_{variants}, S_i, g_1, g_2, \ldots, g_{34}, c_1, c_2, c_3, c_4 \right]
$$

| Feature Set | Dimension | Description |
|-------------|-----------|-------------|
| $n_{variants}$ | 1 | Total pathogenic variant count |
| $S_i$ | 1 | Composite severity score |
| $g_1, \ldots, g_{34}$ | 34 | Binary gene indicators (1 if patient carries variant in gene) |
| $c_1, \ldots, c_4$ | 4 | Consequence counts (missense, nonsense, frameshift, splice) |
| **Total** | **40** | |

### 2.4.3 Federated K-Means Algorithm

The clustering algorithm proceeded iteratively across $P = 5$ population nodes:

**Algorithm: Federated K-Means**

```
Input: K (number of clusters), ε (convergence threshold), T_max (max iterations)
Initialize: Global centroids μ₁⁽⁰⁾, μ₂⁽⁰⁾, ..., μ_K⁽⁰⁾ randomly

For t = 1 to T_max:
    For each population node p = 1 to P:
        # LOCAL STEP 1: Assign patients to nearest centroid
        For each patient i in population p:
            z_i = argmin_k ||x_i - μ_k⁽ᵗ⁻¹⁾||²
        
        # LOCAL STEP 2: Compute local cluster statistics
        For k = 1 to K:
            n_k⁽ᵖ⁾ = |{i : z_i = k}|                    # Local cluster size
            μ_k⁽ᵖ⁾ = (1/n_k⁽ᵖ⁾) Σ_{z_i=k} x_i          # Local centroid
        
        Send {n_k⁽ᵖ⁾, μ_k⁽ᵖ⁾} to central coordinator
    
    # GLOBAL AGGREGATION
    For k = 1 to K:
        μ_k⁽ᵗ⁾ = Σ_p (n_k⁽ᵖ⁾ × μ_k⁽ᵖ⁾) / Σ_p n_k⁽ᵖ⁾    # Weighted average
    
    # CONVERGENCE CHECK
    Δ = Σ_k ||μ_k⁽ᵗ⁾ - μ_k⁽ᵗ⁻¹⁾||²
    If Δ < ε: break

Output: Cluster assignments {z_i} and final centroids {μ_k}
```

**Hyperparameters:**
- $K = 5$ (determined by elbow method and silhouette analysis)
- $\epsilon = 10^{-6}$
- $T_{max} = 100$
- Distance metric: Euclidean

### 2.4.4 Cluster Number Selection

The optimal cluster count was determined using two complementary methods:

**Elbow Method (Within-Cluster Sum of Squares):**

$$
\text{WCSS}(K) = \sum_{k=1}^{K} \sum_{i \in C_k} \|\mathbf{x}_i - \boldsymbol{\mu}_k\|^2
$$

**Silhouette Score:**

For each patient *i*:
$$
s(i) = \frac{b(i) - a(i)}{\max\{a(i), b(i)\}}
$$

Where:
- $a(i)$ = mean distance to other patients in same cluster
- $b(i)$ = mean distance to patients in nearest other cluster

Mean silhouette score across patients indicated optimal separation at $K = 5$.

### 2.4.5 Cluster Interpretation

The resulting five clusters exhibited distinct severity and progression profiles:

| Cluster | n | Mean Severity Score | Modal Progression | % Modal | Interpretation |
|---------|---|--------------------|--------------------|---------|----------------|
| C0 | 287 | 5.00 | Slow | 74.2% (213/287) | Mild |
| C1 | 3,287 | 7.27 | Moderate | 92.9% (3,055/3,287) | Moderate-A |
| C2 | 1,999 | 9.30 | Fast | 96.5% (1,930/1,999) | Severe |
| C3 | 8,957 | 0.00 | N/A | 100% | Control (Healthy) |
| C4 | 845 | 6.32 | Moderate | 91.2% (771/845) | Moderate-B |

**Key Finding:** Clusters C1 and C4, while both classified as "Moderate" severity, emerged as distinct molecular subtypes, providing evidence for **severity-stratified molecular heterogeneity** within the moderate disease category.

**Severity Gradient Validation:**

The cluster severity ordering (C0 < C4 < C1 < C2) was validated by:
1. Monotonic increase in mean severity scores
2. Concordant shift in predicted progression rates
3. Distinct gene enrichment patterns (detailed in Section 2.6)

---

## 2.5 Phase 2: Mechanistic Pathway Annotation

### 2.5.1 Pathway Ontology Definition

Seven canonical ALS-associated molecular pathways were defined based on comprehensive literature review of motor neuron degeneration mechanisms:

| Pathway | Code | Primary Biological Process |
|---------|------|---------------------------|
| Proteostasis | PROT | Protein folding quality control, autophagy, ubiquitin-proteasome system |
| RNA Metabolism | RNA | Pre-mRNA splicing, stress granule dynamics, nucleocytoplasmic transport |
| Cytoskeletal/Axonal Transport | CYTO | Microtubule dynamics, motor protein function, neurofilament organization |
| Mitochondrial Dysfunction | MITO | Oxidative phosphorylation, reactive oxygen species homeostasis, apoptosis |
| Excitotoxicity | EXCITO | Glutamatergic neurotransmission, calcium homeostasis, AMPA/NMDA receptor function |
| Vesicle Trafficking | VES | Endosomal sorting, autophagosome-lysosome fusion, synaptic vesicle cycling |
| DNA Damage Response | DNA | Genome stability, DNA repair, R-loop resolution |

### 2.5.2 Gene-Pathway Mapping

Each of the 34 ALS-associated genes was mapped to one or more pathways based on established molecular mechanisms. The mapping schema permitted **pleiotropy**, reflecting biological reality where single genes participate in multiple cellular processes.

**Complete Gene-Pathway Assignment:**

| Gene | Pathway(s) | Mechanistic Basis |
|------|------------|-------------------|
| *SOD1* | PROT, MITO, EXCITO | Misfolded aggregates disrupt proteostasis; mitochondrial localization impairs respiration; EAAT2 cleavage reduces glutamate clearance |
| *C9orf72* | PROT, RNA, EXCITO | Dipeptide repeat proteins inhibit proteasome; RNA foci sequester splicing factors; AMPA receptor upregulation |
| *TARDBP* | RNA, EXCITO | TDP-43 mislocalization disrupts splicing; ADAR2 downregulation increases Ca²⁺-permeable AMPA receptors |
| *FUS* | RNA, MITO | Nuclear export defects; interaction with ATP synthase impairs mitochondrial function |
| *VCP* | PROT | Autophagosome maturation failure |
| *UBQLN2* | PROT | Ubiquitin-proteasome dysfunction |
| *OPTN* | PROT | Autophagy receptor dysfunction |
| *SQSTM1* | PROT | p62-mediated selective autophagy impairment |
| *TBK1* | PROT | Autophagy initiation kinase deficiency |
| *CCNF* | PROT | E3 ubiquitin ligase dysfunction |
| *MATR3* | RNA | mRNA nuclear export block |
| *HNRNPA1* | RNA | Stress granule dysregulation |
| *HNRNPA2B1* | RNA | Stress granule dysregulation |
| *ANG* | RNA | Ribosomal biogenesis impairment |
| *ELP3* | RNA | tRNA modification defects |
| *TUBA4A* | CYTO | Microtubule destabilization |
| *PFN1* | CYTO | Actin polymerization failure |
| *NEFH* | CYTO | Neurofilament accumulation |
| *DCTN1* | CYTO | Retrograde axonal transport failure |
| *KIF5A* | CYTO | Anterograde axonal transport failure |
| *CHCHD10* | MITO | Cristae disruption, mtDNA instability |
| *SIGMAR1* | MITO | MAM calcium dysregulation |
| *ATXN2* | MITO | NADPH oxidase activation, ROS surge |
| *C19orf12* | MITO | Mitochondrial iron dysregulation |
| *ALS2* | VES | Endosome-lysosome fusion failure |
| *CHMP2B* | VES | ESCRT-III dysfunction |
| *VAPB* | VES | ER-Golgi transport defects |
| *FIG4* | VES | Lysosomal biogenesis failure |
| *SPG11* | VES, DNA | Lysosome reformation failure; DNA repair defects |
| *NEK1* | DNA | DNA damage checkpoint failure |
| *C21orf2* | DNA | DDR defects via NEK1 interaction |
| *SETX* | DNA | R-loop accumulation, genomic instability |
| *UNC13A* | EXCITO | Synaptic vesicle release defects |
| *DAO* | EXCITO | D-serine metabolism, NMDA modulation |

### 2.5.3 Pathway Scoring Algorithm

For each patient *i*, two complementary pathway metrics were computed:

**Binary Pathway Indicator:**

$$
B_{i,p} = \begin{cases}
1 & \text{if } \exists \text{ gene } g \in G_p \text{ such that patient } i \text{ carries variant in } g \\
0 & \text{otherwise}
\end{cases}
$$

Where $G_p$ denotes the set of genes mapped to pathway $p$.

**Pathway Burden Score:**

$$
\text{Score}_{i,p} = \left| \{g \in G_p : V_{i,g} > 0\} \right|
$$

Where $V_{i,g}$ is the count of pathogenic variants in gene $g$ for patient $i$.

This scoring scheme quantifies the **breadth of pathway disruption** by counting unique genes affected, rather than total variant count, thereby avoiding inflation from multiple variants in a single gene.

### 2.5.4 Summary Metrics

**Derived Variables per Patient:**

| Variable | Computation | Interpretation |
|----------|-------------|----------------|
| `n_pathways_affected` | $\sum_{p=1}^{7} B_{i,p}$ | Number of distinct pathways with ≥1 affected gene |
| `primary_pathway` | $\arg\max_p \text{Score}_{i,p}$ | Pathway with highest gene burden |
| `pathway_X_genes` | Concatenated gene list | Genes contributing to pathway X |

**Cohort-Level Summary:**

| Metric | Value |
|--------|-------|
| Total patients | 15,000 |
| Carriers (≥1 pathway affected) | 6,043 (40.3%) |
| Mean pathways per carrier | 2.1 ± 1.3 |
| Patients with ≥3 pathways | 1,847 (30.6% of carriers) |

---

## 2.6 Phase 3: Statistical Analysis of Pathway Patterns

### 2.6.1 Pathway Prevalence by Cluster

Pathway prevalence was computed within each non-control cluster:

$$
\text{Prevalence}_{p,k} = \frac{\sum_{i \in C_k} B_{i,p}}{|C_k|} \times 100\%
$$

Where $C_k$ denotes the set of carrier patients assigned to cluster $k$.

### 2.6.2 Co-occurrence Analysis

#### 2.6.2.1 Co-occurrence Matrix Construction

For each cluster, a symmetric $7 \times 7$ co-occurrence matrix $\mathbf{M}$ was constructed where:

$$
M_{p,q} = \sum_{i=1}^{n} \mathbb{1}\left[B_{i,p} = 1 \land B_{i,q} = 1\right]
$$

Diagonal elements represent within-pathway counts: $M_{p,p} = \sum_{i} B_{i,p}$

#### 2.6.2.2 Co-occurrence Percentage

To normalize for pathway prevalence asymmetry, co-occurrence percentage was computed relative to the smaller pathway:

$$
\text{Co-occ\%}_{p,q} = \frac{M_{p,q}}{\min(M_{p,p}, M_{q,q})} \times 100\%
$$

This formulation ensures that if all patients with the less common pathway also have the more common pathway, Co-occ% = 100%.

#### 2.6.2.3 Odds Ratio Calculation

For each pathway pair $(p, q)$, association strength was quantified using the odds ratio derived from the $2 \times 2$ contingency table:

|  | Pathway $q$ = 1 | Pathway $q$ = 0 | Total |
|--|-----------------|-----------------|-------|
| **Pathway $p$ = 1** | $a$ | $b$ | $a+b$ |
| **Pathway $p$ = 0** | $c$ | $d$ | $c+d$ |
| **Total** | $a+c$ | $b+d$ | $n$ |

**Odds Ratio:**
$$
\text{OR}_{p,q} = \frac{a \times d}{b \times c}
$$

**95% Confidence Interval** (Woolf logit method):
$$
\ln(\text{OR}) \pm 1.96 \times \sqrt{\frac{1}{a} + \frac{1}{b} + \frac{1}{c} + \frac{1}{d}}
$$

$$
\text{CI}_{95\%} = \left[ \exp\left(\ln(\text{OR}) - 1.96 \times \text{SE}\right), \exp\left(\ln(\text{OR}) + 1.96 \times \text{SE}\right) \right]
$$

**Effect Size Interpretation:**

| Odds Ratio | Effect Size Classification |
|------------|---------------------------|
| OR > 3.0 | Large positive association |
| 1.5 < OR ≤ 3.0 | Medium positive association |
| 0.67 ≤ OR ≤ 1.5 | Negligible association |
| 0.33 ≤ OR < 0.67 | Medium negative association |
| OR < 0.33 | Large negative association |

#### 2.6.2.4 Jaccard Similarity Index

To quantify pathway overlap independent of marginal frequencies:

$$
J_{p,q} = \frac{|A_p \cap A_q|}{|A_p \cup A_q|} = \frac{M_{p,q}}{M_{p,p} + M_{q,q} - M_{p,q}}
$$

Where $A_p$ = set of patients with pathway $p$ affected.

Jaccard index ranges from 0 (no overlap) to 1 (complete overlap).

### 2.6.3 Correlation Analysis

#### 2.6.3.1 Spearman Rank Correlation

Spearman's ρ was computed between pathway burden scores to assess monotonic (potentially non-linear) relationships:

$$
\rho_{p,q} = 1 - \frac{6 \sum_{i=1}^{n} d_i^2}{n(n^2 - 1)}
$$

Where $d_i = R(\text{Score}_{i,p}) - R(\text{Score}_{i,q})$ is the difference in ranks for patient $i$.

#### 2.6.3.2 Pearson Correlation

Pearson's $r$ was computed to specifically detect linear dose-response relationships:

$$
r_{p,q} = \frac{\sum_{i=1}^{n}(\text{Score}_{i,p} - \overline{\text{Score}_p})(\text{Score}_{i,q} - \overline{\text{Score}_q})}{\sqrt{\sum_{i=1}^{n}(\text{Score}_{i,p} - \overline{\text{Score}_p})^2} \sqrt{\sum_{i=1}^{n}(\text{Score}_{i,q} - \overline{\text{Score}_q})^2}}
$$

**Rationale for Dual Analysis:**

| Method | Detects | Use Case |
|--------|---------|----------|
| Spearman ρ | Monotonic relationships | Robust to outliers, non-normal distributions |
| Pearson $r$ | Linear relationships | Tests dose-dependent biological hypothesis |

**Linearity Assessment:**

The difference $|\rho - r|$ was computed for each pathway pair. Small differences indicate approximately linear relationships.

*Observed in overall cohort:*
- Maximum $|\rho - r|$ = 0.034
- Mean $|\rho - r|$ = 0.009

The close agreement (Δ < 0.05 for all pairs) confirmed that pathway relationships are **approximately linear**, supporting Pearson correlation for dose-dependent inference.

#### 2.6.3.3 Correlation Strength Classification

| Correlation | Strength | Biological Interpretation |
|-------------|----------|--------------------------|
| $|r| \geq 0.7$ | Strong | Linear dose-dependent relationship |
| $0.3 \leq |r| < 0.7$ | Moderate | Meaningful association with variability |
| $|r| < 0.3$ | Weak | Minimal linear relationship |

### 2.6.4 Quadrant Classification of Pathway Relationships

Pathway pairs were categorized into biological relationship types based on the joint distribution of co-occurrence frequency and correlation strength:

**Classification Criteria:**

| Category | Criteria | Biological Interpretation |
|----------|----------|--------------------------|
| **Dose-dependent** | Co-occ% > 50% AND $r$ > 0.5 | Cascading failure: when one pathway intensifies, the other scales proportionally |
| **Threshold effect** | Co-occ% > 50% AND $|r|$ < 0.3 | Binary disruption: pathways co-occur but intensities are independent |
| **Distinct subtypes** | Co-occ% < 30% AND $r$ < −0.3 | Mutually exclusive: suggests separate molecular subtypes |
| **Other** | Not meeting above criteria | Mixed or transitional patterns |

**Exemplar Classifications (Overall Cohort):**

| Pathway Pair | Co-occ% | Pearson $r$ | Classification |
|--------------|---------|-------------|----------------|
| Mitochondrial × Excitotoxicity | 86.5% | 0.687 | Dose-dependent |
| RNA Metabolism × Mitochondrial | 60.4% | 0.214 | Threshold effect |
| Proteostasis × Vesicle Trafficking | 24.4% | −0.408 | Distinct subtypes |

### 2.6.5 Cross-Cluster Statistical Comparison

#### 2.6.5.1 Chi-Square Test for Prevalence Heterogeneity

For each pathway, a $4 \times 2$ contingency table (4 clusters × pathway present/absent) was analyzed:

$$
\chi^2 = \sum_{i=1}^{4} \sum_{j=1}^{2} \frac{(O_{ij} - E_{ij})^2}{E_{ij}}
$$

Where:
- $O_{ij}$ = observed count in cell $(i,j)$
- $E_{ij} = \frac{R_i \times C_j}{N}$ = expected count under independence
- Degrees of freedom: $df = (4-1)(2-1) = 3$

#### 2.6.5.2 Kruskal-Wallis Test for Score Differences

For pathway burden scores across clusters:

$$
H = \frac{12}{N(N+1)} \sum_{k=1}^{K} \frac{R_k^2}{n_k} - 3(N+1)
$$

Where:
- $N$ = total sample size
- $K$ = number of clusters
- $R_k$ = sum of ranks in cluster $k$
- $n_k$ = sample size of cluster $k$

Under $H_0$, $H$ approximately follows $\chi^2$ distribution with $df = K - 1$.

#### 2.6.5.3 Effect Size Measures

Given the large sample size ($n = 6,043$ carriers), statistical significance was anticipated for most comparisons. **Effect sizes were therefore prioritized** for interpretation of clinical and biological relevance.

**Epsilon-Squared ($\varepsilon^2$)** for Kruskal-Wallis:

$$
\varepsilon^2 = \frac{H}{N - 1}
$$

| $\varepsilon^2$ | Effect Size |
|-----------------|-------------|
| ≥ 0.14 | Large |
| 0.06 – 0.14 | Medium |
| 0.01 – 0.06 | Small |
| < 0.01 | Negligible |

**Cohen's $h$** for pairwise prevalence differences:

$$
h = 2 \arcsin\left(\sqrt{p_1}\right) - 2 \arcsin\left(\sqrt{p_2}\right)
$$

Where $p_1$ and $p_2$ are prevalence proportions in the two clusters being compared.

| $|h|$ | Effect Size |
|-------|-------------|
| ≥ 0.80 | Large |
| 0.50 – 0.80 | Medium |
| 0.20 – 0.50 | Small |
| < 0.20 | Negligible |

### 2.6.6 Multiple Testing Correction

For the correlation matrix (21 unique pairwise comparisons among 7 pathways), Bonferroni correction was applied:

$$
\alpha_{adjusted} = \frac{0.05}{21} = 0.0024
$$

Correlations with $p < 0.0024$ were considered statistically significant.

---

## 2.7 Network Analysis

### 2.7.1 Network Construction

Pathway relationships were represented as an undirected weighted graph $G = (V, E, W)$ where:

- **Vertices** ($V$): Seven molecular pathways
- **Edges** ($E$): Pathway pairs with co-occurrence count > 0
- **Weights** ($W$): Co-occurrence counts

### 2.7.2 Visual Encoding

| Element | Encoding | Formula |
|---------|----------|---------|
| Node size | Pathway prevalence | $r = \frac{\sqrt{n_{pathway}}}{5} + 15$ pixels |
| Node color | Pathway identity | Fixed color palette |
| Edge width | Co-occurrence strength | $w = \ln(\text{count} + 1) \times 2$ |
| Edge color | Effect size | Red (large), Orange (medium), Gray (negligible) |

### 2.7.3 Network Layout

Nodes were arranged in circular layout with angular position:

$$
\theta_p = \frac{2\pi \times \text{index}(p)}{7} - \frac{\pi}{2}
$$

Coordinates:
$$
x_p = x_{center} + R \cos(\theta_p)
$$
$$
y_p = y_{center} + R \sin(\theta_p)
$$

Where $R = 200$ pixels (layout radius).

---

## 2.8 Software Environment

### 2.8.1 Computational Stack

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.10.12 | Primary analysis environment |
| pandas | 2.0.3 | Data manipulation |
| NumPy | 1.24.3 | Numerical computation |
| SciPy | 1.11.1 | Statistical tests |
| scikit-learn | 1.3.0 | K-means clustering |
| D3.js | 7.8.5 | Interactive visualization |

### 2.8.2 Analysis Pipeline

The complete pipeline comprised three executable modules:

1. `phase1_federated_clustering.py` — Federated K-means implementation with cross-population aggregation
2. `phase2_pathway_annotation_all.py` — Gene-pathway mapping and patient scoring
3. `stage4_pathway_analysis_clusters.py` — Statistical analysis, correlation matrices, and network generation

### 2.8.3 Code and Data Availability

All analysis scripts and intermediate datasets are available at [repository URL to be inserted]. The ClinVar source data is publicly accessible through NCBI (https://www.ncbi.nlm.nih.gov/clinvar/).

---

## 2.9 Methodological Considerations

### 2.9.1 Sample Size Justification

With $n = 6,043$ carriers, statistical power exceeded 99% for detecting:
- Small correlations ($r = 0.10$) at $\alpha = 0.05$
- Small prevalence differences (Cohen's $h = 0.20$) between clusters

Consequently, effect sizes rather than $p$-values were emphasized for interpretation.

### 2.9.2 Limitations

1. **Synthetic Data**: Patients were computationally simulated based on variant-level characteristics; clinical phenotype heterogeneity and environmental factors were not modeled.

2. **Pathway Boundary Definitions**: Gene-pathway assignments were based on primary literature and may not capture context-dependent or tissue-specific pathway involvement.

3. **Epistatic Effects**: Complex genetic interactions between variants were not explicitly modeled in the severity score.

4. **Population Structure**: While five superpopulations were represented, within-population stratification (e.g., ancestry-informative markers) was not addressed.

5. **Temporal Dynamics**: The analysis represents a cross-sectional snapshot; longitudinal pathway evolution was not modeled.

---

## 2.10 Ethical Considerations

This study utilized exclusively synthetic data generated from publicly available ClinVar variant annotations. No human subjects were enrolled, and no identifiable patient information was used. Institutional review board approval was not required.

---

*Statistical analyses were performed in Python 3.10 using scipy.stats for inferential tests. Interactive visualizations were developed using D3.js v7.8.5. All analyses were conducted between January and February 2026.*
