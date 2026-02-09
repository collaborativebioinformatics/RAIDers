#!/usr/bin/env python3
"""
Phase 2: Pathway Annotation Script - WEIGHTED SCORING + CLUSTER EMBEDDING
==============================================================
Annotates patients with pathway scores using:
  1. Gene Evidence Weights (penetrance, literature support)
  2. Gene-Gene Interaction Terms (synergy/redundancy)
  3. Normalized Pathway Burden Scores

Author: Lynn Kim
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
"""

import pandas as pd
import numpy as np
from pathlib import Path
import glob
import os
from collections import defaultdict

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FOLDER = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/ALS_Synthetic_Data/phase2_output_weighted")
OUTPUT_FOLDER = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/ALS_Synthetic_Data/phase2_output_weighted_clusters")

# =============================================================================
# PATHWAY DEFINITIONS
# =============================================================================

PATHWAYS = [
    'Proteostasis',
    'RNA_Metabolism',
    'Cytoskeletal_Axonal_Transport',
    'Mitochondrial',
    'Excitotoxicity',
    'Vesicle_Trafficking',
    'DNA_Damage',
]

# Pathway sizes for normalization (number of known genes)
PATHWAY_SIZES = {
    'Proteostasis': 8,
    'RNA_Metabolism': 7,
    'Cytoskeletal_Axonal_Transport': 5,
    'Mitochondrial': 6,
    'Excitotoxicity': 5,
    'Vesicle_Trafficking': 5,
    'DNA_Damage': 4,
}

# =============================================================================
# PHASE 1 CLUSTER DEFINITIONS
# =============================================================================
# Clusters derived from Federated K-means in Phase 1.
# Assignment is deterministic based on severity_score ranges.

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


# =============================================================================
# GENE EVIDENCE WEIGHTS
# =============================================================================
# 
# Weight = Penetrance × Evidence_Multiplier × Pathway_Centrality
#
# Penetrance: 
#   1.0 = High (>80% of carriers develop ALS)
#   0.7 = Moderate (40-80%)
#   0.5 = Variable (20-40%)
#   0.3 = Low/Risk modifier (<20%)
#
# Evidence_Multiplier:
#   1.5 = Strong (>50 publications, clear functional mechanism)
#   1.0 = Moderate (10-50 publications)
#   0.7 = Limited (<10 publications)
#   0.5 = Emerging (few reports, mechanism unclear)
#
# Pathway_Centrality:
#   1.5 = Core/hub gene (central to pathway function)
#   1.0 = Standard member
#   0.7 = Peripheral member
#   0.5 = Indirect/modifier
#

GENE_WEIGHTS = {
    # =========================================================================
    # HIGH CONFIDENCE CAUSATIVE GENES
    # =========================================================================
    'SOD1': {
        'penetrance': 1.0,      # High penetrance, well-established
        'evidence': 1.5,        # >1000 publications, clear mechanism
        'centrality': 1.5,      # Hub gene for proteostasis
        'weight': 2.25,         # 1.0 × 1.5 × 1.5
        'confidence': 'Definitive',
        'notes': 'First ALS gene identified (1993), extensive functional studies'
    },
    'C9ORF72': {
        'penetrance': 0.7,      # Incomplete penetrance, variable expressivity
        'evidence': 1.5,        # >500 publications since 2011
        'centrality': 1.5,      # Affects multiple pathways
        'weight': 1.575,        # 0.7 × 1.5 × 1.5
        'confidence': 'Definitive',
        'notes': 'Most common genetic cause of fALS/sALS in Europeans'
    },
    'TARDBP': {
        'penetrance': 1.0,      # High penetrance
        'evidence': 1.5,        # Extensive TDP-43 literature
        'centrality': 1.5,      # Central to RNA metabolism
        'weight': 2.25,         # 1.0 × 1.5 × 1.5
        'confidence': 'Definitive',
        'notes': 'TDP-43 pathology present in >95% of ALS cases'
    },
    'FUS': {
        'penetrance': 1.0,      # High penetrance, especially juvenile ALS
        'evidence': 1.5,        # Well-characterized
        'centrality': 1.5,      # Hub for RNA processing
        'weight': 2.25,         # 1.0 × 1.5 × 1.5
        'confidence': 'Definitive',
        'notes': 'Associated with aggressive juvenile-onset ALS'
    },
    
    # =========================================================================
    # MODERATE CONFIDENCE CAUSATIVE GENES
    # =========================================================================
    'VCP': {
        'penetrance': 0.7,      # Variable, often with IBM or FTD
        'evidence': 1.0,        # Good evidence
        'centrality': 1.0,      # Important for proteostasis
        'weight': 0.7,          # 0.7 × 1.0 × 1.0
        'confidence': 'Strong',
        'notes': 'Multisystem proteinopathy, ALS in ~10% of carriers'
    },
    'TBK1': {
        'penetrance': 0.5,      # Haploinsufficiency, variable
        'evidence': 1.0,        # Growing evidence
        'centrality': 1.2,      # Key kinase for autophagy
        'weight': 0.6,          # 0.5 × 1.0 × 1.2
        'confidence': 'Strong',
        'notes': 'Activates OPTN and p62, critical for autophagy'
    },
    'OPTN': {
        'penetrance': 0.5,      # Variable penetrance
        'evidence': 1.0,        
        'centrality': 1.0,      
        'weight': 0.5,          # 0.5 × 1.0 × 1.0
        'confidence': 'Strong',
        'notes': 'Autophagy receptor, TBK1 substrate'
    },
    'SQSTM1': {
        'penetrance': 0.5,      
        'evidence': 1.0,        
        'centrality': 1.0,      
        'weight': 0.5,          # 0.5 × 1.0 × 1.0
        'confidence': 'Strong',
        'notes': 'p62, autophagy receptor'
    },
    'UBQLN2': {
        'penetrance': 0.7,      # X-linked, higher in males
        'evidence': 1.0,        
        'centrality': 1.0,      
        'weight': 0.7,          # 0.7 × 1.0 × 1.0
        'confidence': 'Strong',
        'notes': 'X-linked, may show anticipation'
    },
    'CCNF': {
        'penetrance': 0.5,      
        'evidence': 0.7,        # Limited studies
        'centrality': 0.7,      
        'weight': 0.245,        # 0.5 × 0.7 × 0.7
        'confidence': 'Moderate',
        'notes': 'E3 ubiquitin ligase complex'
    },
    
    # =========================================================================
    # RNA METABOLISM GENES
    # =========================================================================
    'MATR3': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.35,         # 0.5 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'Nuclear matrix protein, mRNA export'
    },
    'HNRNPA1': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 0.7,      
        'weight': 0.245,        # 0.5 × 0.7 × 0.7
        'confidence': 'Moderate',
        'notes': 'Stress granule component'
    },
    'HNRNPA2B1': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 0.7,      
        'weight': 0.245,        # 0.5 × 0.7 × 0.7
        'confidence': 'Moderate',
        'notes': 'Stress granule component, similar to HNRNPA1'
    },
    'ANG': {
        'penetrance': 0.3,      # Low penetrance risk factor
        'evidence': 1.0,        
        'centrality': 0.7,      
        'weight': 0.21,         # 0.3 × 1.0 × 0.7
        'confidence': 'Moderate',
        'notes': 'Angiogenin, ribosomal biogenesis'
    },
    'ELP3': {
        'penetrance': 0.3,      
        'evidence': 0.5,        
        'centrality': 0.7,      
        'weight': 0.105,        # 0.3 × 0.5 × 0.7
        'confidence': 'Limited',
        'notes': 'tRNA modification, elongator complex'
    },
    
    # =========================================================================
    # CYTOSKELETAL / AXONAL TRANSPORT GENES
    # =========================================================================
    'TUBA4A': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.35,         # 0.5 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'Tubulin, microtubule network'
    },
    'PFN1': {
        'penetrance': 0.7,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.49,         # 0.7 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'Profilin, actin dynamics'
    },
    'NEFH': {
        'penetrance': 0.3,      # Risk modifier
        'evidence': 1.0,        
        'centrality': 0.7,      
        'weight': 0.21,         # 0.3 × 1.0 × 0.7
        'confidence': 'Moderate',
        'notes': 'Neurofilament heavy chain'
    },
    'DCTN1': {
        'penetrance': 0.3,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.21,         # 0.3 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'Dynactin, retrograde transport'
    },
    'KIF5A': {
        'penetrance': 0.5,      
        'evidence': 1.0,        
        'centrality': 1.0,      
        'weight': 0.5,          # 0.5 × 1.0 × 1.0
        'confidence': 'Strong',
        'notes': 'Kinesin, anterograde transport'
    },
    
    # =========================================================================
    # MITOCHONDRIAL GENES
    # =========================================================================
    'CHCHD10': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.35,         # 0.5 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'Mitochondrial cristae morphology'
    },
    'SIGMAR1': {
        'penetrance': 0.3,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.21,         # 0.3 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'MAM Ca2+ signaling'
    },
    'ATXN2': {
        'penetrance': 0.3,      # Risk modifier (intermediate repeats)
        'evidence': 1.5,        # Well-characterized as modifier
        'centrality': 0.7,      
        'weight': 0.315,        # 0.3 × 1.5 × 0.7
        'confidence': 'Strong (as modifier)',
        'notes': 'PolyQ intermediate repeats increase ALS risk'
    },
    'C19orf12': {
        'penetrance': 0.3,      
        'evidence': 0.5,        
        'centrality': 0.7,      
        'weight': 0.105,        # 0.3 × 0.5 × 0.7
        'confidence': 'Limited',
        'notes': 'Mitochondrial iron metabolism'
    },
    
    # =========================================================================
    # VESICLE TRAFFICKING GENES
    # =========================================================================
    'ALS2': {
        'penetrance': 0.7,      # Juvenile ALS, recessive
        'evidence': 1.0,        
        'centrality': 1.0,      
        'weight': 0.7,          # 0.7 × 1.0 × 1.0
        'confidence': 'Strong',
        'notes': 'Alsin, Rab5 GEF, juvenile ALS'
    },
    'CHMP2B': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.35,         # 0.5 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'ESCRT-III, endosomal sorting'
    },
    'VAPB': {
        'penetrance': 0.7,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.49,         # 0.7 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'ER-Golgi transport'
    },
    'FIG4': {
        'penetrance': 0.3,      
        'evidence': 0.7,        
        'centrality': 0.7,      
        'weight': 0.147,        # 0.3 × 0.7 × 0.7
        'confidence': 'Limited',
        'notes': 'PI(3,5)P2 metabolism'
    },
    'SPG11': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 0.7,      
        'weight': 0.245,        # 0.5 × 0.7 × 0.7
        'confidence': 'Moderate',
        'notes': 'Spatacsin, lysosome reformation'
    },
    
    # =========================================================================
    # DNA DAMAGE RESPONSE GENES
    # =========================================================================
    'NEK1': {
        'penetrance': 0.3,      # Risk factor, not fully penetrant
        'evidence': 1.0,        # Well-characterized
        'centrality': 1.2,      # Key DDR kinase
        'weight': 0.36,         # 0.3 × 1.0 × 1.2
        'confidence': 'Strong (as risk)',
        'notes': 'DNA damage response, cilia'
    },
    'C21orf2': {
        'penetrance': 0.3,      
        'evidence': 0.5,        
        'centrality': 0.7,      
        'weight': 0.105,        # 0.3 × 0.5 × 0.7
        'confidence': 'Limited',
        'notes': 'NEK1 interactor'
    },
    'SETX': {
        'penetrance': 0.5,      
        'evidence': 0.7,        
        'centrality': 1.0,      
        'weight': 0.35,         # 0.5 × 0.7 × 1.0
        'confidence': 'Moderate',
        'notes': 'Senataxin, R-loop resolution'
    },
    
    # =========================================================================
    # EXCITOTOXICITY GENES
    # =========================================================================
    'UNC13A': {
        'penetrance': 0.3,      # Risk modifier
        'evidence': 1.0,        
        'centrality': 0.7,      
        'weight': 0.21,         # 0.3 × 1.0 × 0.7
        'confidence': 'Moderate',
        'notes': 'Synaptic vesicle release'
    },
    'DAO': {
        'penetrance': 0.3,      
        'evidence': 0.5,        
        'centrality': 0.7,      
        'weight': 0.105,        # 0.3 × 0.5 × 0.7
        'confidence': 'Limited',
        'notes': 'D-amino acid oxidase, NMDA modulation'
    },
}

# Default weight for genes not in our curated list
DEFAULT_GENE_WEIGHT = {
    'penetrance': 0.3,
    'evidence': 0.5,
    'centrality': 0.5,
    'weight': 0.075,  # 0.3 × 0.5 × 0.5
    'confidence': 'Unknown',
    'notes': 'Gene not in curated ALS gene list'
}

# =============================================================================
# GENE-GENE INTERACTION MATRIX
# =============================================================================
#
# Synergy multiplier > 1.0: Genes amplify each other's effect
# Redundancy multiplier < 1.0: Genes have overlapping function
# Neutral = 1.0: Independent effects
#
# Format: (gene1, gene2): multiplier
#

GENE_INTERACTIONS = {
    # =========================================================================
    # SYNERGISTIC INTERACTIONS (multiplier > 1.0)
    # =========================================================================
    
    # TDP-43 and FUS: Both RNA-binding proteins that aggregate, compete for targets
    ('TARDBP', 'FUS'): {
        'multiplier': 1.3,
        'type': 'synergy',
        'mechanism': 'Both aggregate in cytoplasm, sequester overlapping RNA targets',
        'evidence': 'Multiple studies show co-aggregation and synthetic lethality in models'
    },
    
    # OPTN and TBK1: TBK1 phosphorylates OPTN to activate autophagy
    ('OPTN', 'TBK1'): {
        'multiplier': 1.5,
        'type': 'synergy',
        'mechanism': 'TBK1 activates OPTN; dual loss = complete autophagy failure',
        'evidence': 'Biochemical and genetic studies confirm epistatic relationship'
    },
    
    # SOD1 and SIGMAR1: Both affect mitochondrial Ca2+ handling
    ('SOD1', 'SIGMAR1'): {
        'multiplier': 1.4,
        'type': 'synergy',
        'mechanism': 'SOD1 aggregates in mitochondria; SIGMAR1 regulates MAM Ca2+',
        'evidence': 'Both converge on mitochondrial dysfunction, amplified toxicity'
    },
    
    # C9ORF72 and TARDBP: DPR proteins enhance TDP-43 aggregation
    ('C9ORF72', 'TARDBP'): {
        'multiplier': 1.3,
        'type': 'synergy',
        'mechanism': 'C9 DPR proteins impair nuclear transport, enhancing TDP-43 mislocalization',
        'evidence': 'Studies show C9 expansion accelerates TDP-43 pathology'
    },
    
    # C9ORF72 and FUS: Similar mechanism as above
    ('C9ORF72', 'FUS'): {
        'multiplier': 1.2,
        'type': 'synergy',
        'mechanism': 'DPR proteins affect nuclear transport of FUS',
        'evidence': 'Moderate evidence for interaction'
    },
    
    # SOD1 and mitochondrial genes
    ('SOD1', 'CHCHD10'): {
        'multiplier': 1.3,
        'type': 'synergy',
        'mechanism': 'Both impair mitochondrial function through different mechanisms',
        'evidence': 'Convergent mitochondrial toxicity'
    },
    
    # FUS and ATP synthase (CHCHD10 related)
    ('FUS', 'CHCHD10'): {
        'multiplier': 1.2,
        'type': 'synergy',
        'mechanism': 'FUS disrupts ATP synthase; CHCHD10 affects cristae',
        'evidence': 'Both target mitochondrial energy production'
    },
    
    # DNA damage and proteostasis
    ('NEK1', 'TBK1'): {
        'multiplier': 1.2,
        'type': 'synergy',
        'mechanism': 'Both are kinases that regulate stress responses',
        'evidence': 'Limited but suggestive evidence'
    },
    
    # =========================================================================
    # REDUNDANT INTERACTIONS (multiplier < 1.0)
    # =========================================================================
    
    # HNRNPA1 and HNRNPA2B1: Same stress granule mechanism
    ('HNRNPA1', 'HNRNPA2B1'): {
        'multiplier': 0.7,
        'type': 'redundancy',
        'mechanism': 'Both have prion-like domains, same stress granule role',
        'evidence': 'Highly similar function, partial redundancy expected'
    },
    
    # OPTN and SQSTM1: Both autophagy receptors
    ('OPTN', 'SQSTM1'): {
        'multiplier': 0.8,
        'type': 'redundancy',
        'mechanism': 'Both are autophagy receptors with overlapping cargo',
        'evidence': 'Functional redundancy in autophagy receptor function'
    },
    
    # DCTN1 and KIF5A: Opposite transport directions
    ('DCTN1', 'KIF5A'): {
        'multiplier': 0.9,
        'type': 'partial_redundancy',
        'mechanism': 'Retrograde vs anterograde transport, some functional overlap',
        'evidence': 'Both required for axonal transport but different directions'
    },
    
    # TBK1 and SQSTM1: TBK1 activates p62
    ('TBK1', 'SQSTM1'): {
        'multiplier': 0.85,
        'type': 'epistatic',
        'mechanism': 'TBK1 phosphorylates p62; mutations may be in same pathway',
        'evidence': 'Epistatic relationship, not fully additive'
    },
}

# =============================================================================
# CROSS-PATHWAY AMPLIFICATION
# =============================================================================
# When a gene affects multiple pathways, its damage creates feedback loops

CROSS_PATHWAY_BONUS = {
    # Gene: bonus multiplier for each additional pathway affected
    'SOD1': 1.15,      # Strong cross-pathway effects
    'C9ORF72': 1.15,   # Affects proteostasis, RNA, excitotoxicity
    'TARDBP': 1.1,     # RNA + excitotoxicity
    'FUS': 1.1,        # RNA + mitochondrial
    'SPG11': 1.05,     # Vesicle + DNA damage
}

DEFAULT_CROSS_PATHWAY_BONUS = 1.05  # Small bonus for any multi-pathway gene

# =============================================================================
# GENE-PATHWAY MAPPING (with mechanisms)
# =============================================================================

GENE_PATHWAY_MECHANISM = {
    'SOD1': {
        'pathways': ['Proteostasis', 'Mitochondrial', 'Excitotoxicity'],
        'mechanism_short': 'Protein misfolding; mitochondrial toxicity; glutamate transporter cleavage',
    },
    'C9ORF72': {
        'pathways': ['Proteostasis', 'RNA_Metabolism', 'Excitotoxicity'],
        'mechanism_short': 'DPR toxicity; RNA foci; AMPA receptor dysregulation',
    },
    'VCP': {
        'pathways': ['Proteostasis'],
        'mechanism_short': 'Autophagosome maturation failure; TDP-43 accumulation',
    },
    'UBQLN2': {
        'pathways': ['Proteostasis'],
        'mechanism_short': 'Ubiquitin-proteasome clearance failure',
    },
    'OPTN': {
        'pathways': ['Proteostasis'],
        'mechanism_short': 'Autophagy receptor dysfunction',
    },
    'SQSTM1': {
        'pathways': ['Proteostasis'],
        'mechanism_short': 'Autophagy receptor dysfunction (p62)',
    },
    'TBK1': {
        'pathways': ['Proteostasis'],
        'mechanism_short': 'Autophagy initiation kinase deficiency',
    },
    'CCNF': {
        'pathways': ['Proteostasis'],
        'mechanism_short': 'E3 ubiquitin-ligase dysfunction',
    },
    'TARDBP': {
        'pathways': ['RNA_Metabolism', 'Excitotoxicity'],
        'mechanism_short': 'Nuclear-cytoplasmic mislocalization; ADAR2 downregulation',
    },
    'FUS': {
        'pathways': ['RNA_Metabolism', 'Mitochondrial'],
        'mechanism_short': 'Nuclear-cytoplasmic mislocalization; ATP synthase disruption',
    },
    'MATR3': {
        'pathways': ['RNA_Metabolism'],
        'mechanism_short': 'mRNA nuclear export block',
    },
    'HNRNPA1': {
        'pathways': ['RNA_Metabolism'],
        'mechanism_short': 'Stress granule dysregulation',
    },
    'HNRNPA2B1': {
        'pathways': ['RNA_Metabolism'],
        'mechanism_short': 'Stress granule dysregulation',
    },
    'ANG': {
        'pathways': ['RNA_Metabolism'],
        'mechanism_short': 'Ribosomal biogenesis impairment',
    },
    'ELP3': {
        'pathways': ['RNA_Metabolism'],
        'mechanism_short': 'tRNA modification defects',
    },
    'TUBA4A': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism_short': 'Microtubule destabilization',
    },
    'PFN1': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism_short': 'Actin polymerization failure',
    },
    'NEFH': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism_short': 'Neurofilament accumulation',
    },
    'DCTN1': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism_short': 'Retrograde transport failure',
    },
    'KIF5A': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism_short': 'Anterograde transport failure',
    },
    'CHCHD10': {
        'pathways': ['Mitochondrial'],
        'mechanism_short': 'Cristae disruption; mtDNA instability',
    },
    'SIGMAR1': {
        'pathways': ['Mitochondrial'],
        'mechanism_short': 'MAM calcium dysregulation',
    },
    'ATXN2': {
        'pathways': ['Mitochondrial'],
        'mechanism_short': 'NADPH oxidase activation; ROS surge',
    },
    'C19orf12': {
        'pathways': ['Mitochondrial'],
        'mechanism_short': 'Mitochondrial iron dysregulation',
    },
    'ALS2': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism_short': 'Endosome-lysosome fusion failure',
    },
    'CHMP2B': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism_short': 'ESCRT-III dysfunction',
    },
    'VAPB': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism_short': 'ER-Golgi transport defects',
    },
    'FIG4': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism_short': 'Lysosomal biogenesis failure',
    },
    'SPG11': {
        'pathways': ['Vesicle_Trafficking', 'DNA_Damage'],
        'mechanism_short': 'Lysosome reformation failure; DNA repair defects',
    },
    'NEK1': {
        'pathways': ['DNA_Damage'],
        'mechanism_short': 'DNA damage checkpoint failure',
    },
    'C21orf2': {
        'pathways': ['DNA_Damage'],
        'mechanism_short': 'NEK1 interaction; DDR defects',
    },
    'SETX': {
        'pathways': ['DNA_Damage'],
        'mechanism_short': 'R-loop accumulation; genomic instability',
    },
    'UNC13A': {
        'pathways': ['Excitotoxicity'],
        'mechanism_short': 'Synaptic transmission defects',
    },
    'DAO': {
        'pathways': ['Excitotoxicity'],
        'mechanism_short': 'D-serine metabolism; NMDA modulation',
    },
}

# =============================================================================
# SCORING FUNCTIONS
# =============================================================================

def parse_genes(all_genes_str):
    """Parse the 'all_genes' column into a list of gene names."""
    if pd.isna(all_genes_str) or all_genes_str == '':
        return []
    if ';' in str(all_genes_str):
        genes = [g.strip() for g in str(all_genes_str).split(';')]
    else:
        genes = [g.strip() for g in str(all_genes_str).split(',')]
    return [g.upper() for g in genes if g]


def get_gene_weight(gene_name):
    """Get the evidence weight for a gene."""
    gene_clean = gene_name.upper().strip()
    if gene_clean in GENE_WEIGHTS:
        return GENE_WEIGHTS[gene_clean]['weight']
    return DEFAULT_GENE_WEIGHT['weight']


def get_gene_info(gene_name):
    """Get full gene information including weight and confidence."""
    gene_clean = gene_name.upper().strip()
    if gene_clean in GENE_WEIGHTS:
        return GENE_WEIGHTS[gene_clean]
    return DEFAULT_GENE_WEIGHT


def calculate_interaction_multiplier(genes_in_pathway):
    """
    Calculate the interaction multiplier for a set of genes in a pathway.
    
    For each pair of genes, check if there's a known interaction.
    Multiply all interaction multipliers together.
    """
    if len(genes_in_pathway) < 2:
        return 1.0
    
    multiplier = 1.0
    gene_list = list(genes_in_pathway)
    
    for i in range(len(gene_list)):
        for j in range(i + 1, len(gene_list)):
            g1, g2 = sorted([gene_list[i], gene_list[j]])
            pair = (g1, g2)
            
            if pair in GENE_INTERACTIONS:
                multiplier *= GENE_INTERACTIONS[pair]['multiplier']
    
    return multiplier


def calculate_cross_pathway_bonus(gene, n_pathways):
    """
    Calculate bonus for genes affecting multiple pathways.
    
    Genes that affect multiple pathways may create feedback loops.
    """
    if n_pathways <= 1:
        return 1.0
    
    gene_clean = gene.upper().strip()
    bonus_per_pathway = CROSS_PATHWAY_BONUS.get(gene_clean, DEFAULT_CROSS_PATHWAY_BONUS)
    
    # Apply bonus for each additional pathway beyond the first
    return bonus_per_pathway ** (n_pathways - 1)


def get_patient_pathway_scores(patient_genes):
    """
    Calculate comprehensive pathway scores for a single patient.
    
    Returns dict with multiple scoring approaches:
    1. Binary (0/1)
    2. Simple count
    3. Weighted sum
    4. With interaction multiplier
    5. Normalized burden score
    """
    unique_genes = set(patient_genes)
    
    # Track genes and weights per pathway
    pathway_data = {pathway: {
        'genes': set(),
        'weights': [],
        'mechanisms': []
    } for pathway in PATHWAYS}
    
    # Map genes to pathways
    for gene in unique_genes:
        gene_clean = gene.upper().strip()
        gene_info = GENE_PATHWAY_MECHANISM.get(gene_clean)
        
        if gene_info:
            weight = get_gene_weight(gene_clean)
            n_pathways = len(gene_info['pathways'])
            cross_bonus = calculate_cross_pathway_bonus(gene_clean, n_pathways)
            
            for pathway in gene_info['pathways']:
                pathway_data[pathway]['genes'].add(gene_clean)
                # Apply cross-pathway bonus to the weight
                pathway_data[pathway]['weights'].append(weight * cross_bonus)
                pathway_data[pathway]['mechanisms'].append(
                    f"{gene_clean}: {gene_info['mechanism_short']}"
                )
    
    result = {}
    total_burden = 0
    
    for pathway in PATHWAYS:
        data = pathway_data[pathway]
        genes_in_pathway = data['genes']
        n_genes = len(genes_in_pathway)
        
        # Binary
        result[f'pathway_{pathway}'] = 1 if n_genes > 0 else 0
        
        # Simple count (original method)
        result[f'pathway_{pathway}_count'] = n_genes
        
        # Weighted sum
        weighted_sum = sum(data['weights'])
        result[f'pathway_{pathway}_weighted'] = round(weighted_sum, 3)
        
        # Interaction multiplier
        interaction_mult = calculate_interaction_multiplier(genes_in_pathway)
        result[f'pathway_{pathway}_interaction'] = round(interaction_mult, 3)
        
        # Normalized burden score
        # Formula: (weighted_sum × interaction_multiplier) / sqrt(pathway_size)
        pathway_size = PATHWAY_SIZES[pathway]
        if weighted_sum > 0:
            burden = (weighted_sum * interaction_mult) / np.sqrt(pathway_size)
        else:
            burden = 0
        result[f'pathway_{pathway}_burden'] = round(burden, 3)
        total_burden += burden
        
        # Gene list and mechanisms
        result[f'pathway_{pathway}_genes'] = '; '.join(sorted(genes_in_pathway)) if genes_in_pathway else ''
        result[f'pathway_{pathway}_mechanisms'] = ' | '.join(data['mechanisms']) if data['mechanisms'] else ''
    
    # Summary columns
    n_pathways = sum(1 for p in PATHWAYS if result[f'pathway_{p}'] == 1)
    result['n_pathways_affected'] = n_pathways
    
    # Find primary pathway (highest burden)
    if n_pathways > 0:
        primary = max(PATHWAYS, key=lambda p: result[f'pathway_{p}_burden'])
        result['primary_pathway'] = primary
    else:
        result['primary_pathway'] = 'None'
    
    # Total and composite scores
    result['total_burden_score'] = round(total_burden, 3)
    
    # Multi-pathway adjustment
    multi_factor = 1 + 0.1 * (n_pathways - 1) if n_pathways > 0 else 1
    result['multi_pathway_factor'] = round(multi_factor, 2)
    result['composite_score'] = round(total_burden * multi_factor, 3)
    
    # Gene weight details for transparency
    gene_weights_detail = []
    for gene in sorted(unique_genes):
        gene_clean = gene.upper().strip()
        info = get_gene_info(gene_clean)
        gene_weights_detail.append(f"{gene_clean}={info['weight']:.3f}")
    result['gene_weights_detail'] = '; '.join(gene_weights_detail) if gene_weights_detail else ''
    
    return result


def calculate_population_percentiles(df, score_column='composite_score'):
    """
    Calculate percentile scores relative to the population.
    
    Carriers only - percentile indicates where patient falls
    relative to other carriers.
    """
    carriers = df[df['n_pathways_affected'] > 0].copy()
    
    if len(carriers) == 0:
        df['burden_percentile'] = np.nan
        return df
    
    # Calculate percentile for each carrier
    df['burden_percentile'] = df[score_column].apply(
        lambda x: round(100 * (carriers[score_column] < x).sum() / len(carriers), 1)
        if x > 0 else 0
    )
    
    return df


def assign_risk_category(composite_score):
    """
    Assign genetic risk category based on composite score.
    
    Thresholds based on expected distribution:
    - Low: < 1.0 (single low-penetrance variant)
    - Moderate: 1.0 - 2.5 (multiple low variants or one moderate)
    - High: 2.5 - 5.0 (high penetrance gene or multiple moderates)
    - Very High: > 5.0 (multiple high penetrance genes or strong interactions)
    """
    if composite_score < 0.5:
        return 'Minimal'
    elif composite_score < 1.0:
        return 'Low'
    elif composite_score < 2.5:
        return 'Moderate'
    elif composite_score < 5.0:
        return 'High'
    else:
        return 'Very High'


# =============================================================================
# FILE PROCESSING
# =============================================================================


# =============================================================================
# PHASE 1 CLUSTER ASSIGNMENT FUNCTIONS
# =============================================================================

def assign_phase1_cluster(severity_score):
    """Assign Phase 1 cluster based on severity_score ranges."""
    if severity_score == 0:
        return 3  # Unaffected
    elif 4.0 <= severity_score <= 5.7:
        return 0  # Mild
    elif 6.0 <= severity_score <= 6.6:
        return 4  # Moderate-ALS2
    elif 6.9 <= severity_score <= 7.8:
        return 1  # Moderate-SOD1
    elif 8.4 <= severity_score <= 10.0:
        return 2  # Severe
    else:
        return -1  # Unmatched


def add_cluster_columns(df):
    """Add Phase 1 cluster and cluster_label columns to dataframe."""
    df = df.copy()
    df['cluster'] = df['severity_score'].apply(assign_phase1_cluster)
    df['cluster_label'] = df['cluster'].map(CLUSTER_LABELS)
    
    # Reorder columns to place cluster info after severity_category
    cols = df.columns.tolist()
    if 'severity_category' in cols and 'cluster' in cols:
        insert_pos = cols.index('severity_category') + 1
        cols.remove('cluster')
        cols.remove('cluster_label')
        cols.insert(insert_pos, 'cluster')
        cols.insert(insert_pos + 1, 'cluster_label')
        df = df[cols]
    return df


def print_cluster_summary(df, title="PHASE 1 CLUSTER SUMMARY"):
    """Print cluster distribution and validation."""
    print(f"\n--- {title} ---")
    print("\nCluster Distribution:")
    for cid in sorted(df['cluster'].unique()):
        n = len(df[df['cluster'] == cid])
        label = CLUSTER_LABELS.get(cid, 'Unknown')
        pct = 100 * n / len(df)
        print(f"  {cid} ({label}): {n} ({pct:.1f}%)")
    
    print("\nCluster vs Severity Category:")
    print(pd.crosstab(df['cluster'], df['severity_category']))


def process_population_file(filepath):
    """Process a single population CSV file."""
    
    print(f"\nProcessing: {filepath}")
    df = pd.read_csv(filepath)
    print(f"  Loaded {len(df)} patients")
    
    # Extract population code
    filename = os.path.basename(filepath)
    pop_code = filename.replace('client_', '').replace('.csv', '')
    df['superpopulation'] = pop_code
    
    # Columns to keep
    base_cols = [
        'patient_id', 'superpopulation', 'n_variants_carried', 'primary_gene',
        'severity_score', 'severity_category', 'predicted_progression', 'all_genes'
    ]
    keep_cols = [c for c in base_cols if c in df.columns]
    
    # Calculate pathway scores
    print("  Calculating weighted pathway scores...")
    pathway_data = []
    
    for idx, row in df.iterrows():
        genes = parse_genes(row.get('all_genes', ''))
        scores = get_patient_pathway_scores(genes)
        pathway_data.append(scores)
    
    pathway_df = pd.DataFrame(pathway_data)
    result_df = pd.concat([df[keep_cols].reset_index(drop=True), pathway_df], axis=1)
    
    # Add risk category
    result_df['genetic_risk_category'] = result_df['composite_score'].apply(assign_risk_category)
    
    # Calculate percentiles
    result_df = calculate_population_percentiles(result_df)
    
    carriers = result_df[result_df['n_pathways_affected'] > 0]
    print(f"  Carriers: {len(carriers)}")
    print(f"  Mean composite score (carriers): {carriers['composite_score'].mean():.3f}")
    
    # Add Phase 1 cluster columns
    result_df = add_cluster_columns(result_df)
    
    return result_df


def main():
    print("=" * 70)
    print("Phase 2: Pathway Annotation - WEIGHTED SCORING VERSION")
    print("=" * 70)
    
    # Print methodology summary
    print("""
SCORING METHODOLOGY:
  1. Gene Evidence Weight = Penetrance × Evidence × Centrality
  2. Interaction Multiplier = Product of pairwise gene interactions
  3. Pathway Burden = (Σ Weights × Interaction) / √(pathway_size)
  4. Composite Score = Total Burden × Multi-pathway Factor
    """)
    
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    print(f"\nInput folder:  {INPUT_FOLDER}")
    print(f"Output folder: {OUTPUT_FOLDER}")
    
    # Find input files
    pattern = os.path.join(INPUT_FOLDER, "client_*.csv")
    all_files = glob.glob(pattern)
    input_files = [f for f in all_files if '_carriers' not in f and '_pathways' not in f]
    
    print(f"\nFound {len(input_files)} population files:")
    for f in sorted(input_files):
        print(f"  - {os.path.basename(f)}")
    
    if len(input_files) == 0:
        print("\nERROR: No client_*.csv files found!")
        return
    
    # Process each population
    all_dfs = []
    
    for filepath in sorted(input_files):
        df = process_population_file(filepath)
        all_dfs.append(df)
        
        pop_code = os.path.basename(filepath).replace('client_', '').replace('.csv', '')
        output_path = os.path.join(OUTPUT_FOLDER, f"client_{pop_code}_weighted.csv")
        df.to_csv(output_path, index=False)
        print(f"  Saved: {output_path}")
    
    # Combine all populations
    print("\n" + "=" * 70)
    print("Combining all populations...")
    combined_df = pd.concat(all_dfs, ignore_index=True)
    
    # Recalculate percentiles on combined data
    combined_df = calculate_population_percentiles(combined_df)
    
    # Validate cluster assignments
    print("\n" + "=" * 70)
    print("PHASE 1 CLUSTER EMBEDDING VALIDATION")
    print("=" * 70)
    print_cluster_summary(combined_df)
    
    unassigned = combined_df[combined_df['cluster'] == -1]
    if len(unassigned) > 0:
        print(f"\n⚠️ WARNING: {len(unassigned)} patients with unmatched severity_score!")
    else:
        print("\n✓ All patients assigned to Phase 1 clusters!")
    
    combined_path = os.path.join(OUTPUT_FOLDER, "patients_with_pathways_weighted_clusters.csv")
    combined_df.to_csv(combined_path, index=False)
    print(f"\nSaved: {combined_path}")
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    carriers = combined_df[combined_df['n_pathways_affected'] > 0]
    
    print(f"\nTotal patients: {len(combined_df)}")
    print(f"Total carriers: {len(carriers)}")
    
    print("\n--- Composite Score Distribution (Carriers) ---")
    print(f"  Mean:   {carriers['composite_score'].mean():.3f}")
    print(f"  Median: {carriers['composite_score'].median():.3f}")
    print(f"  Std:    {carriers['composite_score'].std():.3f}")
    print(f"  Min:    {carriers['composite_score'].min():.3f}")
    print(f"  Max:    {carriers['composite_score'].max():.3f}")
    
    print("\n--- Risk Category Distribution ---")
    for cat in ['Minimal', 'Low', 'Moderate', 'High', 'Very High']:
        count = len(carriers[carriers['genetic_risk_category'] == cat])
        pct = 100 * count / len(carriers) if len(carriers) > 0 else 0
        print(f"  {cat}: {count} ({pct:.1f}%)")
    
    print("\n--- Pathway Burden Scores (Mean, Carriers) ---")
    for pathway in PATHWAYS:
        col = f'pathway_{pathway}_burden'
        affected = carriers[carriers[f'pathway_{pathway}'] == 1]
        if len(affected) > 0:
            mean_burden = affected[col].mean()
            print(f"  {pathway}: {mean_burden:.3f} (n={len(affected)})")
    
    # Save methodology documentation
    methodology_path = os.path.join(OUTPUT_FOLDER, "SCORING_METHODOLOGY.md")
    with open(methodology_path, 'w') as f:
        f.write(__doc__)
        f.write("\n\n## Gene Weights Table\n\n")
        f.write("| Gene | Penetrance | Evidence | Centrality | Weight | Confidence |\n")
        f.write("|------|------------|----------|------------|--------|------------|\n")
        for gene, info in sorted(GENE_WEIGHTS.items()):
            f.write(f"| {gene} | {info['penetrance']} | {info['evidence']} | {info['centrality']} | {info['weight']:.3f} | {info['confidence']} |\n")
        
        f.write("\n\n## Gene-Gene Interactions\n\n")
        f.write("| Gene 1 | Gene 2 | Multiplier | Type | Mechanism |\n")
        f.write("|--------|--------|------------|------|------------|\n")
        for (g1, g2), info in sorted(GENE_INTERACTIONS.items()):
            f.write(f"| {g1} | {g2} | {info['multiplier']} | {info['type']} | {info['mechanism'][:50]}... |\n")
    
    print(f"\nSaved methodology documentation: {methodology_path}")
    
    print("\n" + "=" * 70)
    print("COMPARISON: OLD vs NEW SCORING")
    print("=" * 70)
    
    # Show example patients to illustrate difference
    sample = carriers.nlargest(5, 'composite_score')[
        ['patient_id', 'all_genes', 'n_pathways_affected', 
         'pathway_Proteostasis_count', 'pathway_Proteostasis_burden',
         'composite_score', 'genetic_risk_category']
    ]
    print("\nTop 5 patients by composite score:")
    print(sample.to_string())
    
    # Pathway distribution by cluster
    print("\n" + "=" * 70)
    print("PATHWAY DISTRIBUTION BY CLUSTER")
    print("=" * 70)
    
    for cid in [0, 1, 2, 4]:
        cdf = carriers[carriers['cluster'] == cid]
        n = len(cdf)
        if n == 0:
            continue
        label = CLUSTER_LABELS[cid]
        print(f"\n{label} (n={n}):")
        for pathway in PATHWAYS:
            col = f'pathway_{pathway}'
            if col in cdf.columns:
                count = cdf[col].sum()
                pct = 100 * count / n
                print(f"  {pathway}: {count} ({pct:.1f}%)")
    
    print("\n✓ Phase 2 weighted annotation with cluster embedding complete!")


if __name__ == "__main__":
    main()
