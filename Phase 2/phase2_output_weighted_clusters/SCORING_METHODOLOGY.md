
Phase 2: Pathway Annotation Script - WEIGHTED SCORING + CLUSTER EMBEDDING
==============================================================
Annotates patients with pathway scores using:
  1. Gene Evidence Weights (penetrance, literature support)
  2. Gene-Gene Interaction Terms (synergy/redundancy)
  3. Normalized Pathway Burden Scores

Author: [Your name]
Version: 2.0 (Weighted)

================================================================================
METHODOLOGY DOCUMENTATION
================================================================================

This script implements three layers of scoring to improve upon simple gene counts:

┌─────────────────────────────────────────────────────────────────────────────┐
│ LAYER 1: GENE EVIDENCE WEIGHT (W_gene)                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│ Each gene is assigned a weight based on ALS literature evidence:            │
│                                                                             │
│   W_gene = Penetrance × Evidence_Multiplier × Pathway_Centrality            │
│                                                                             │
│ Where:                                                                      │
│   • Penetrance (0.1 - 1.0): Probability mutation causes disease             │
│     - High penetrance (>80%): SOD1, FUS, TARDBP = 1.0                       │
│     - Moderate (40-80%): C9ORF72, VCP = 0.7                                 │
│     - Low (<40%): Risk variants like ATXN2, UNC13A = 0.3                    │
│                                                                             │
│   • Evidence_Multiplier (0.5 - 1.5): Literature support strength            │
│     - Strong (>50 publications, functional studies): 1.5                    │
│     - Moderate (10-50 publications): 1.0                                    │
│     - Limited (<10 publications): 0.5                                       │
│                                                                             │
│   • Pathway_Centrality (0.5 - 1.5): How central to pathway function         │
│     - Core/hub gene: 1.5                                                    │
│     - Standard member: 1.0                                                  │
│     - Peripheral: 0.5                                                       │
│                                                                             │
│ FORMULA: W_gene = penetrance × evidence × centrality                        │
│                                                                             │
│ Example: SOD1 = 1.0 × 1.5 × 1.5 = 2.25                                      │
│          ATXN2 = 0.3 × 1.0 × 1.0 = 0.30                                     │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ LAYER 2: GENE-GENE INTERACTION SCORE (I_interaction)                        │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│ When multiple genes affect the same pathway, they may interact:             │
│                                                                             │
│ A) SYNERGISTIC INTERACTIONS (multiplier > 1.0)                              │
│    Genes that amplify each other's effects                                  │
│                                                                             │
│    Known synergies in ALS:                                                  │
│    • TDP-43 + FUS: Both aggregate, compete for RNA targets (1.3×)           │
│    • SOD1 + SIGMAR1: Both impair mitochondrial Ca2+ (1.4×)                  │
│    • C9ORF72 + TARDBP: DPR proteins enhance TDP-43 aggregation (1.3×)       │
│    • OPTN + TBK1: TBK1 phosphorylates OPTN; both lost = autophagy           │
│      collapse (1.5×)                                                        │
│                                                                             │
│ B) REDUNDANT INTERACTIONS (multiplier < 1.0)                                │
│    Genes in same sub-pathway = diminishing returns                          │
│                                                                             │
│    • HNRNPA1 + HNRNPA2B1: Same stress granule mechanism (0.7×)              │
│    • OPTN + SQSTM1: Both autophagy receptors, partly redundant (0.8×)       │
│    • DCTN1 + KIF5A: Opposite transport directions, less overlap (0.9×)      │
│                                                                             │
│ C) CROSS-PATHWAY AMPLIFICATION                                              │
│    Genes affecting multiple pathways may create feedback loops:             │
│                                                                             │
│    • SOD1 (Proteostasis + Mitochondrial): Misfolded protein in              │
│      mitochondria creates ROS → more misfolding (1.2× cross-pathway)        │
│    • FUS (RNA_Metabolism + Mitochondrial): RNA dysregulation affects        │
│      mitochondrial-encoded proteins (1.2× cross-pathway)                    │
│                                                                             │
│ FORMULA:                                                                    │
│   I_interaction = Π (interaction_multiplier for each gene pair)             │
│                                                                             │
│ Example: Patient with OPTN + TBK1 + SOD1 in Proteostasis                    │
│   I = 1.5 (OPTN×TBK1) × 1.0 (no SOD1 interactions) = 1.5                    │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ LAYER 3: PATHWAY BURDEN SCORE (B_pathway)                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│ Raw weighted scores are normalized to account for pathway size:             │
│                                                                             │
│ A) RAW WEIGHTED SCORE:                                                      │
│    S_raw = Σ (W_gene for each gene in pathway) × I_interaction              │
│                                                                             │
│ B) NORMALIZED BURDEN SCORE:                                                 │
│    B_pathway = S_raw / sqrt(N_genes_in_pathway)                             │
│                                                                             │
│    Using sqrt normalization because:                                        │
│    - Linear (÷N) over-penalizes large pathways                              │
│    - No normalization ignores pathway size entirely                         │
│    - Sqrt is standard in gene set enrichment (GSEA-style)                   │
│                                                                             │
│ C) PERCENTILE SCORE (0-100):                                                │
│    For clinical interpretation, burden scores are converted to              │
│    percentiles based on the population distribution                         │
│                                                                             │
│ PATHWAY SIZES (for reference):                                              │
│    Proteostasis: 8 genes → sqrt = 2.83                                      │
│    RNA_Metabolism: 7 genes → sqrt = 2.65                                    │
│    Cytoskeletal_Axonal_Transport: 5 genes → sqrt = 2.24                     │
│    Mitochondrial: 6 genes (including cross-pathway) → sqrt = 2.45           │
│    Excitotoxicity: 5 genes (including cross-pathway) → sqrt = 2.24          │
│    Vesicle_Trafficking: 5 genes → sqrt = 2.24                               │
│    DNA_Damage: 4 genes → sqrt = 2.00                                        │
│                                                                             │
│ FORMULA:                                                                    │
│   B_pathway = [Σ(W_gene) × I_interaction] / sqrt(pathway_size)              │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ LAYER 4: COMPOSITE PATIENT SCORE                                            │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│ A) TOTAL BURDEN SCORE:                                                      │
│    Total = Σ (B_pathway for all affected pathways)                          │
│                                                                             │
│ B) MULTI-PATHWAY PENALTY/BONUS:                                             │
│    Patients with multiple pathways affected may have worse prognosis        │
│                                                                             │
│    Multi_factor = 1 + 0.1 × (n_pathways - 1)                                │
│                                                                             │
│    Example: 3 pathways → 1 + 0.1 × 2 = 1.2× multiplier                      │
│                                                                             │
│ C) FINAL COMPOSITE SCORE:                                                   │
│    Composite = Total × Multi_factor                                         │
│                                                                             │
│ INTERPRETATION:                                                             │
│    0-1.0: Low genetic burden                                                │
│    1.0-2.5: Moderate genetic burden                                         │
│    2.5-5.0: High genetic burden                                             │
│    >5.0: Very high genetic burden                                           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘

================================================================================
OUTPUT COLUMNS
================================================================================

For each pathway:
  • pathway_{name}                 : Binary (0/1)
  • pathway_{name}_count           : Simple gene count (original method)
  • pathway_{name}_weighted        : Sum of gene weights
  • pathway_{name}_interaction     : Interaction multiplier
  • pathway_{name}_burden          : Normalized burden score
  • pathway_{name}_genes           : Gene names
  • pathway_{name}_mechanisms      : Mechanism descriptions

Summary columns:
  • total_burden_score            : Sum of all pathway burdens
  • composite_score               : With multi-pathway adjustment
  • burden_percentile             : Population percentile (0-100)
  • genetic_risk_category         : Low/Moderate/High/Very High

================================================================================


## Gene Weights Table

| Gene | Penetrance | Evidence | Centrality | Weight | Confidence |
|------|------------|----------|------------|--------|------------|
| ALS2 | 0.7 | 1.0 | 1.0 | 0.700 | Strong |
| ANG | 0.3 | 1.0 | 0.7 | 0.210 | Moderate |
| ATXN2 | 0.3 | 1.5 | 0.7 | 0.315 | Strong (as modifier) |
| C19orf12 | 0.3 | 0.5 | 0.7 | 0.105 | Limited |
| C21orf2 | 0.3 | 0.5 | 0.7 | 0.105 | Limited |
| C9ORF72 | 0.7 | 1.5 | 1.5 | 1.575 | Definitive |
| CCNF | 0.5 | 0.7 | 0.7 | 0.245 | Moderate |
| CHCHD10 | 0.5 | 0.7 | 1.0 | 0.350 | Moderate |
| CHMP2B | 0.5 | 0.7 | 1.0 | 0.350 | Moderate |
| DAO | 0.3 | 0.5 | 0.7 | 0.105 | Limited |
| DCTN1 | 0.3 | 0.7 | 1.0 | 0.210 | Moderate |
| ELP3 | 0.3 | 0.5 | 0.7 | 0.105 | Limited |
| FIG4 | 0.3 | 0.7 | 0.7 | 0.147 | Limited |
| FUS | 1.0 | 1.5 | 1.5 | 2.250 | Definitive |
| HNRNPA1 | 0.5 | 0.7 | 0.7 | 0.245 | Moderate |
| HNRNPA2B1 | 0.5 | 0.7 | 0.7 | 0.245 | Moderate |
| KIF5A | 0.5 | 1.0 | 1.0 | 0.500 | Strong |
| MATR3 | 0.5 | 0.7 | 1.0 | 0.350 | Moderate |
| NEFH | 0.3 | 1.0 | 0.7 | 0.210 | Moderate |
| NEK1 | 0.3 | 1.0 | 1.2 | 0.360 | Strong (as risk) |
| OPTN | 0.5 | 1.0 | 1.0 | 0.500 | Strong |
| PFN1 | 0.7 | 0.7 | 1.0 | 0.490 | Moderate |
| SETX | 0.5 | 0.7 | 1.0 | 0.350 | Moderate |
| SIGMAR1 | 0.3 | 0.7 | 1.0 | 0.210 | Moderate |
| SOD1 | 1.0 | 1.5 | 1.5 | 2.250 | Definitive |
| SPG11 | 0.5 | 0.7 | 0.7 | 0.245 | Moderate |
| SQSTM1 | 0.5 | 1.0 | 1.0 | 0.500 | Strong |
| TARDBP | 1.0 | 1.5 | 1.5 | 2.250 | Definitive |
| TBK1 | 0.5 | 1.0 | 1.2 | 0.600 | Strong |
| TUBA4A | 0.5 | 0.7 | 1.0 | 0.350 | Moderate |
| UBQLN2 | 0.7 | 1.0 | 1.0 | 0.700 | Strong |
| UNC13A | 0.3 | 1.0 | 0.7 | 0.210 | Moderate |
| VAPB | 0.7 | 0.7 | 1.0 | 0.490 | Moderate |
| VCP | 0.7 | 1.0 | 1.0 | 0.700 | Strong |


## Gene-Gene Interactions

| Gene 1 | Gene 2 | Multiplier | Type | Mechanism |
|--------|--------|------------|------|------------|
| C9ORF72 | FUS | 1.2 | synergy | DPR proteins affect nuclear transport of FUS... |
| C9ORF72 | TARDBP | 1.3 | synergy | C9 DPR proteins impair nuclear transport, enhancin... |
| DCTN1 | KIF5A | 0.9 | partial_redundancy | Retrograde vs anterograde transport, some function... |
| FUS | CHCHD10 | 1.2 | synergy | FUS disrupts ATP synthase; CHCHD10 affects cristae... |
| HNRNPA1 | HNRNPA2B1 | 0.7 | redundancy | Both have prion-like domains, same stress granule ... |
| NEK1 | TBK1 | 1.2 | synergy | Both are kinases that regulate stress responses... |
| OPTN | SQSTM1 | 0.8 | redundancy | Both are autophagy receptors with overlapping carg... |
| OPTN | TBK1 | 1.5 | synergy | TBK1 activates OPTN; dual loss = complete autophag... |
| SOD1 | CHCHD10 | 1.3 | synergy | Both impair mitochondrial function through differe... |
| SOD1 | SIGMAR1 | 1.4 | synergy | SOD1 aggregates in mitochondria; SIGMAR1 regulates... |
| TARDBP | FUS | 1.3 | synergy | Both aggregate in cytoplasm, sequester overlapping... |
| TBK1 | SQSTM1 | 0.85 | epistatic | TBK1 phosphorylates p62; mutations may be in same ... |
