# Weighted Pathway Scoring Methodology

## Overview

This document describes the weighted scoring system for ALS pathway analysis, replacing simple gene counts with evidence-based weights, interaction terms, and normalized burden scores.

---

## Table of Contents

1. [Why Weighted Scoring?](#1-why-weighted-scoring)
2. [Layer 1: Gene Evidence Weights](#2-layer-1-gene-evidence-weights)
3. [Layer 2: Gene-Gene Interactions](#3-layer-2-gene-gene-interactions)
4. [Layer 3: Pathway Burden Scores](#4-layer-3-pathway-burden-scores)
5. [Layer 4: Composite Patient Score](#5-layer-4-composite-patient-score)
6. [Complete Calculation Example](#6-complete-calculation-example)
7. [Output Columns Reference](#7-output-columns-reference)
8. [Limitations and Caveats](#8-limitations-and-caveats)

---

## 1. Why Weighted Scoring?

### Problems with Simple Gene Counts

| Issue | Example | Impact |
|-------|---------|--------|
| No penetrance weighting | SOD1 (>90% penetrance) counts same as ATXN2 (risk modifier) | Overestimates risk for low-penetrance variants |
| No interaction modeling | OPTN + TBK1 together = autophagy collapse, but count = 2 | Misses synergistic effects |
| No pathway normalization | Proteostasis (8 genes) vs DNA_Damage (4 genes) | Larger pathways always score higher |
| Ignores evidence strength | Well-studied SOD1 = recently discovered gene | May overweight uncertain associations |

### Weighted Scoring Solution

```
Simple Count:     Patient A (SOD1) = 1,  Patient B (ATXN2) = 1  → Same score!
Weighted Score:   Patient A (SOD1) = 2.25, Patient B (ATXN2) = 0.32  → Reflects biology
```

---

## 2. Layer 1: Gene Evidence Weights

### Formula

```
W_gene = Penetrance × Evidence_Multiplier × Pathway_Centrality
```

### Component Definitions

#### Penetrance (0.1 - 1.0)
Probability that a mutation carrier develops ALS.

| Value | Category | Definition | Examples |
|-------|----------|------------|----------|
| 1.0 | High | >80% of carriers develop ALS | SOD1, FUS, TARDBP |
| 0.7 | Moderate | 40-80% penetrance | C9ORF72, VCP, UBQLN2 |
| 0.5 | Variable | 20-40% penetrance | OPTN, TBK1, KIF5A |
| 0.3 | Low/Modifier | <20%, risk modifier | ATXN2, UNC13A, NEK1 |

**Source**: ClinGen ALS gene curation, literature meta-analyses

#### Evidence Multiplier (0.5 - 1.5)
Strength of literature support for the gene-ALS association.

| Value | Category | Definition | Examples |
|-------|----------|------------|----------|
| 1.5 | Strong | >50 publications, clear mechanism, functional studies | SOD1, C9ORF72, TDP-43 |
| 1.0 | Moderate | 10-50 publications, established association | VCP, OPTN, KIF5A |
| 0.7 | Limited | <10 publications, mechanism partially understood | CCNF, MATR3 |
| 0.5 | Emerging | Few reports, mechanism unclear | ELP3, C21orf2 |

**Source**: PubMed publication counts, ALSoD database, ClinGen reviews

#### Pathway Centrality (0.5 - 1.5)
How central the gene is to pathway function.

| Value | Category | Definition | Examples |
|-------|----------|------------|----------|
| 1.5 | Hub/Core | Essential for pathway function, many interactions | SOD1 (proteostasis), TARDBP (RNA) |
| 1.2 | Important | Key regulator, moderate connectivity | TBK1, NEK1 |
| 1.0 | Standard | Typical pathway member | OPTN, CHCHD10 |
| 0.7 | Peripheral | Accessory function, fewer interactions | ANG, FIG4 |
| 0.5 | Indirect | Modifier effect, not directly in pathway | Some risk variants |

**Source**: Protein-protein interaction networks, pathway databases (KEGG, Reactome)

### Example Calculations

```
SOD1:   1.0 × 1.5 × 1.5 = 2.25  (High penetrance, strong evidence, hub gene)
C9ORF72: 0.7 × 1.5 × 1.5 = 1.575 (Moderate penetrance, strong evidence, hub)
ATXN2:  0.3 × 1.5 × 0.7 = 0.315 (Low penetrance modifier, strong evidence, peripheral)
ELP3:   0.3 × 0.5 × 0.7 = 0.105 (Low penetrance, limited evidence, peripheral)
```

### Complete Gene Weight Table

| Gene | Penetrance | Evidence | Centrality | **Weight** | Confidence |
|------|------------|----------|------------|------------|------------|
| SOD1 | 1.0 | 1.5 | 1.5 | **2.250** | Definitive |
| TARDBP | 1.0 | 1.5 | 1.5 | **2.250** | Definitive |
| FUS | 1.0 | 1.5 | 1.5 | **2.250** | Definitive |
| C9ORF72 | 0.7 | 1.5 | 1.5 | **1.575** | Definitive |
| ALS2 | 0.7 | 1.0 | 1.0 | **0.700** | Strong |
| VCP | 0.7 | 1.0 | 1.0 | **0.700** | Strong |
| UBQLN2 | 0.7 | 1.0 | 1.0 | **0.700** | Strong |
| TBK1 | 0.5 | 1.0 | 1.2 | **0.600** | Strong |
| KIF5A | 0.5 | 1.0 | 1.0 | **0.500** | Strong |
| OPTN | 0.5 | 1.0 | 1.0 | **0.500** | Strong |
| SQSTM1 | 0.5 | 1.0 | 1.0 | **0.500** | Strong |
| PFN1 | 0.7 | 0.7 | 1.0 | **0.490** | Moderate |
| VAPB | 0.7 | 0.7 | 1.0 | **0.490** | Moderate |
| NEK1 | 0.3 | 1.0 | 1.2 | **0.360** | Strong (risk) |
| TUBA4A | 0.5 | 0.7 | 1.0 | **0.350** | Moderate |
| CHCHD10 | 0.5 | 0.7 | 1.0 | **0.350** | Moderate |
| CHMP2B | 0.5 | 0.7 | 1.0 | **0.350** | Moderate |
| MATR3 | 0.5 | 0.7 | 1.0 | **0.350** | Moderate |
| SETX | 0.5 | 0.7 | 1.0 | **0.350** | Moderate |
| ATXN2 | 0.3 | 1.5 | 0.7 | **0.315** | Strong (modifier) |
| CCNF | 0.5 | 0.7 | 0.7 | **0.245** | Moderate |
| HNRNPA1 | 0.5 | 0.7 | 0.7 | **0.245** | Moderate |
| HNRNPA2B1 | 0.5 | 0.7 | 0.7 | **0.245** | Moderate |
| SPG11 | 0.5 | 0.7 | 0.7 | **0.245** | Moderate |
| ANG | 0.3 | 1.0 | 0.7 | **0.210** | Moderate |
| NEFH | 0.3 | 1.0 | 0.7 | **0.210** | Moderate |
| DCTN1 | 0.3 | 0.7 | 1.0 | **0.210** | Moderate |
| UNC13A | 0.3 | 1.0 | 0.7 | **0.210** | Moderate |
| SIGMAR1 | 0.3 | 0.7 | 1.0 | **0.210** | Moderate |
| FIG4 | 0.3 | 0.7 | 0.7 | **0.147** | Limited |
| ELP3 | 0.3 | 0.5 | 0.7 | **0.105** | Limited |
| C19orf12 | 0.3 | 0.5 | 0.7 | **0.105** | Limited |
| C21orf2 | 0.3 | 0.5 | 0.7 | **0.105** | Limited |
| DAO | 0.3 | 0.5 | 0.7 | **0.105** | Limited |
| Unknown | 0.3 | 0.5 | 0.5 | **0.075** | Unknown |

---

## 3. Layer 2: Gene-Gene Interactions

### Concept

When multiple genes affect the same pathway, they may:
- **Synergize** (amplify each other's effect) → multiplier > 1.0
- **Be redundant** (overlapping function) → multiplier < 1.0
- **Act independently** → multiplier = 1.0

### Formula

```
I_interaction = ∏ (interaction_multiplier for each gene pair)
```

For a patient with genes A, B, C in a pathway:
```
I = multiplier(A,B) × multiplier(A,C) × multiplier(B,C)
```

### Synergistic Interactions (multiplier > 1.0)

| Gene 1 | Gene 2 | Multiplier | Mechanism |
|--------|--------|------------|-----------|
| OPTN | TBK1 | **1.50** | TBK1 phosphorylates OPTN; both lost = autophagy collapse |
| SOD1 | SIGMAR1 | **1.40** | Both impair mitochondrial Ca²⁺ handling |
| TARDBP | FUS | **1.30** | Both aggregate, sequester overlapping RNA targets |
| C9ORF72 | TARDBP | **1.30** | DPR proteins enhance TDP-43 mislocalization |
| SOD1 | CHCHD10 | **1.30** | Convergent mitochondrial toxicity |
| C9ORF72 | FUS | **1.20** | DPR proteins affect FUS nuclear transport |
| FUS | CHCHD10 | **1.20** | Both target mitochondrial energy production |
| NEK1 | TBK1 | **1.20** | Both kinases regulate stress responses |

### Redundant Interactions (multiplier < 1.0)

| Gene 1 | Gene 2 | Multiplier | Mechanism |
|--------|--------|------------|-----------|
| HNRNPA1 | HNRNPA2B1 | **0.70** | Same stress granule mechanism, highly similar |
| OPTN | SQSTM1 | **0.80** | Both autophagy receptors, overlapping cargo |
| TBK1 | SQSTM1 | **0.85** | TBK1 activates p62, epistatic relationship |
| DCTN1 | KIF5A | **0.90** | Opposite transport directions, some overlap |

### Example Calculation

Patient with OPTN + TBK1 + SOD1 in Proteostasis pathway:

```
Pairs to check:
  (OPTN, TBK1) = 1.50  ← Known synergy
  (OPTN, SOD1) = 1.00  ← No known interaction
  (TBK1, SOD1) = 1.00  ← No known interaction

I_interaction = 1.50 × 1.00 × 1.00 = 1.50
```

---

## 4. Layer 3: Pathway Burden Scores

### The Problem with Raw Sums

Larger pathways (more known genes) will always have higher raw scores, even if the patient has the same proportion of genes affected.

### Normalization Formula

```
B_pathway = (Σ W_gene × I_interaction) / √(N_genes_in_pathway)
```

**Why square root?**
- Linear normalization (÷N) over-penalizes large pathways
- No normalization ignores pathway size entirely
- Square root is standard in gene set enrichment analysis (GSEA)

### Pathway Sizes

| Pathway | N Genes | √N |
|---------|---------|-----|
| Proteostasis | 8 | 2.83 |
| RNA_Metabolism | 7 | 2.65 |
| Mitochondrial | 6 | 2.45 |
| Cytoskeletal_Axonal_Transport | 5 | 2.24 |
| Excitotoxicity | 5 | 2.24 |
| Vesicle_Trafficking | 5 | 2.24 |
| DNA_Damage | 4 | 2.00 |

### Example Calculation

Patient with SOD1 + TBK1 in Proteostasis:

```
Weighted sum = W(SOD1) + W(TBK1) = 2.25 + 0.60 = 2.85
Interaction  = 1.00 (no known SOD1-TBK1 interaction)
Pathway size = 8 genes, √8 = 2.83

B_proteostasis = (2.85 × 1.00) / 2.83 = 1.007
```

---

## 5. Layer 4: Composite Patient Score

### Total Burden Score

```
Total = Σ (B_pathway for all affected pathways)
```

### Multi-Pathway Adjustment

Patients with multiple pathways affected may have compounding effects:

```
Multi_factor = 1 + 0.1 × (n_pathways - 1)
```

| Pathways | Factor |
|----------|--------|
| 1 | 1.00 |
| 2 | 1.10 |
| 3 | 1.20 |
| 4 | 1.30 |
| 5 | 1.40 |

### Final Composite Score

```
Composite = Total_Burden × Multi_factor
```

### Risk Categories

| Score Range | Category | Interpretation |
|-------------|----------|----------------|
| 0 - 0.5 | Minimal | Likely non-pathogenic or very low risk |
| 0.5 - 1.0 | Low | Single low-penetrance variant |
| 1.0 - 2.5 | Moderate | Multiple low variants OR one moderate |
| 2.5 - 5.0 | High | High penetrance gene OR multiple moderates |
| > 5.0 | Very High | Multiple high penetrance genes OR strong synergies |

---

## 6. Complete Calculation Example

### Patient Profile
- **Genes**: SOD1, SIGMAR1, OPTN
- **Affected Pathways**: Proteostasis, Mitochondrial, Excitotoxicity

### Step 1: Gene Weights

```
W(SOD1)    = 2.250 (high penetrance, definitive)
W(SIGMAR1) = 0.210 (low penetrance, moderate evidence)
W(OPTN)    = 0.500 (variable penetrance, strong evidence)
```

### Step 2: Cross-Pathway Bonus

SOD1 affects 3 pathways → bonus = 1.15^(3-1) = 1.32

```
W(SOD1) adjusted = 2.250 × 1.32 = 2.97
```

### Step 3: Pathway-Specific Calculations

**Proteostasis** (SOD1, OPTN):
```
Weighted sum = 2.97 + 0.50 = 3.47
Interaction  = 1.00 (no known SOD1-OPTN interaction)
Burden       = 3.47 / √8 = 1.227
```

**Mitochondrial** (SOD1, SIGMAR1):
```
Weighted sum = 2.97 + 0.21 = 3.18
Interaction  = 1.40 (SOD1-SIGMAR1 synergy!)
Burden       = (3.18 × 1.40) / √6 = 1.817
```

**Excitotoxicity** (SOD1 only):
```
Weighted sum = 2.97
Interaction  = 1.00 (single gene)
Burden       = 2.97 / √5 = 1.328
```

### Step 4: Composite Score

```
Total Burden = 1.227 + 1.817 + 1.328 = 4.372
Multi-factor = 1 + 0.1 × (3-1) = 1.20
Composite    = 4.372 × 1.20 = 5.246

→ Risk Category: Very High
```

### Comparison to Simple Count

| Method | Score | Interpretation |
|--------|-------|----------------|
| Simple count | 3 genes | No risk stratification |
| Weighted | 5.25 | Very High risk (top percentile) |

---

## 7. Output Columns Reference

### Per-Pathway Columns

| Column | Description | Example |
|--------|-------------|---------|
| `pathway_{name}` | Binary (0/1) | 1 |
| `pathway_{name}_count` | Simple gene count | 2 |
| `pathway_{name}_weighted` | Sum of gene weights | 2.85 |
| `pathway_{name}_interaction` | Interaction multiplier | 1.50 |
| `pathway_{name}_burden` | Normalized burden score | 1.51 |
| `pathway_{name}_genes` | Gene names | "OPTN; TBK1" |
| `pathway_{name}_mechanisms` | Mechanism descriptions | "OPTN: Autophagy receptor..." |

### Summary Columns

| Column | Description | Example |
|--------|-------------|---------|W
| `n_pathways_affected` | Count of pathways with burden > 0 | 3 |
| `primary_pathway` | Pathway with highest burden | "Mitochondrial" |
| `total_burden_score` | Sum of all pathway burdens | 4.37 |
| `multi_pathway_factor` | 1 + 0.1 × (n_pathways - 1) | 1.20 |
| `composite_score` | Total × multi-factor | 5.25 |
| `burden_percentile` | Population percentile (0-100) | 95.2 |
| `genetic_risk_category` | Minimal/Low/Moderate/High/Very High | "Very High" |
| `gene_weights_detail` | Individual gene weights | "SOD1=2.25; OPTN=0.50" |

---

## 8. Limitations and Caveats

### Known Limitations

1. **Weights are estimates**: Penetrance and evidence values are approximations from literature, not precise measurements.

2. **Interaction data is incomplete**: Only well-characterized interactions are included. Unknown interactions default to neutral (1.0).

3. **No variant-level scoring**: All variants in a gene are weighted equally. A severe truncating mutation counts the same as a mild missense.

4. **Population-specific effects not modeled**: Gene frequencies and penetrance may vary by ancestry.

5. **Environmental factors ignored**: Gene-environment interactions not captured.

### Recommended Improvements (Future Work)

1. **Integrate CADD/REVEL scores**: Weight variants by predicted pathogenicity
2. **Add variant-level data**: Distinguish LOF vs missense vs splice variants
3. **Population-specific weights**: Adjust for ancestry-specific penetrance
4. **Machine learning refinement**: Train weights on clinical outcome data
5. **Expand interaction network**: Include protein-protein interaction databases

### Appropriate Use

✅ **Do use for:**
- Research prioritization
- Patient stratification in clinical trials
- Hypothesis generation
- Comparing relative genetic burden

❌ **Do NOT use for:**
- Clinical diagnosis (not validated)
- Predictive genetic testing
- Individual prognosis without clinical correlation

---

## References

1. ClinGen ALS Gene Curation: https://clinicalgenome.org/
2. ALSoD Database: https://alsod.ac.uk/
3. Brown & Al-Chalabi (2017). "ALS Genetics." Nature Reviews Neurology.
4. Taylor et al. (2016). "Deciphering the ALS Molecular Nexus." Science.
5. GSEA Methodology: Subramanian et al. (2005). PNAS.

---

*Document Version: 2.0*  
*Last Updated: January 2025*
