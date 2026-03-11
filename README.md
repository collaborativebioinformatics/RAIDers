# RAIDers (Rare Disease & AI)

<p align="">
  <img width="740" alt="3A821E49-3C43-4990-9FFC-724236CB40C5"
       src="https://github.com/user-attachments/assets/c5f2f452-cb87-4b0c-9bd5-a818cb80408b" />
  
![Python](https://img.shields.io/badge/Python-3.9+-blue)
![Federated Learning](https://img.shields.io/badge/Federated%20Learning-NVFLARE-green)
![Status](https://img.shields.io/badge/Status-Research%20Prototype-yellow)
![License](https://img.shields.io/badge/License-MIT-black)

</p>


**RAIDers** is a federated computational framework for resolving the phenotypic heterogeneity of Amyotrophic Lateral Sclerosis (ALS) while preserving data sovereignty. By integrating global genomic annotations with mechanistic pathway analysis, this framework establishes a scalable architecture for rare disease molecular subtyping.

<img width="6518" height="4170" alt="workflow_flowchart_v3" src="https://github.com/user-attachments/assets/c6cff9cb-5d96-4892-a901-a229b7210376" />

---

## 🎯 Scientific Objectives

- **Subtype Discovery:** Identify coherent molecular signatures across diverse ancestral backgrounds using federated machine learning
- **Pathway-Level Analysis:** Map genetic variants to biological mechanisms to understand disease heterogeneity
- **Federated Feasibility:** Demonstrate that analytical fidelity is maintained when data remains physically separated across institutional nodes
- **Therapeutic Stratification:** Enable precision medicine approaches by identifying distinct molecular subtypes requiring different therapeutic strategies

---

## 📊 Project Phases

| Phase | Status | Description |
|:-----:|:------:|:------------|
| **Phase 1** | ✅ Complete | Federated K-means clustering for molecular subtype discovery |
| **Phase 2** | ✅ Complete | Mechanistic pathway annotation and statistical analysis |
| **Phase 3** | 🔄 In Progress | Multi-institutional collaboration and real-world validation |

---

## 🧬 Phase 1: Federated Subtype Discovery

### Synthetic Cohort Generation

To overcome the "mathematical invisibility" of rare variants, RAIDers employs a digital mutagenesis strategy generating a balanced cohort of **15,000 synthetic patients** partitioned across five ancestral superpopulations (AFR, AMR, EAS, EUR, SAS; n=3,000 each).

**Data Sources:**
- **ClinVar:** ~450 pathogenic ALS variants across 34 genes
- **gnomAD:** Population-specific allele frequencies for ancestral modifiers

### Federated K-Means Clustering

Molecular subtypes are discovered through a privacy-preserving decentralized algorithm:

1. **Local Computation:** Each population node computes cluster assignments using local patient data
2. **Secure Aggregation:** Only centroids (not patient data) are shared with a central coordinator
3. **Global Update:** Federated averaging produces updated global centroids
4. **Convergence:** Process repeats until centroid change < 0.001

### Phase 1 Results

The algorithm identified **5 distinct clusters** with clear clinical stratification:

| Cluster | n | Mean Severity | Progression | Interpretation |
|:-------:|:-:|:-------------:|:-----------:|:---------------|
| C0 | 287 | 5.00 | 74% Slow | Mild |
| C1 | 3,287 | 7.27 | 93% Moderate | Moderate-A (SOD1-like) |
| C2 | 1,999 | 9.30 | 97% Fast | Severe |
| C3 | 8,957 | 0.00 | N/A | Control (No variants) |
| C4 | 845 | 6.32 | 91% Moderate | Moderate-B (ALS2-like) |

**Key Finding:** Clusters 1 and 4 both represent moderate severity but exhibit distinct molecular profiles, providing evidence for **severity-stratified molecular subtypes**.

---

## 🔬 Phase 2: Mechanistic Pathway Analysis

### Pathway Annotation Framework

Phase 2 maps the 34 ALS-associated genes to **7 canonical molecular pathways**:

| Pathway | Genes | Representative |
|:--------|:-----:|:---------------|
| Proteostasis | 8 | *SOD1, C9orf72, VCP, TBK1* |
| RNA Metabolism | 7 | *TARDBP, FUS, MATR3, HNRNPA1* |
| Vesicle Trafficking | 5 | *ALS2, CHMP2B, VAPB, FIG4* |
| Cytoskeletal/Axonal | 5 | *TUBA4A, PFN1, DCTN1, KIF5A* |
| Mitochondrial | 4 | *CHCHD10, SIGMAR1* |
| DNA Damage | 4 | *NEK1, SETX, SPG11* |
| Excitotoxicity | 4 | *UNC13A, DAO* |

### Statistical Analysis Suite

Phase 2 implements comprehensive statistical methods:

- **Co-occurrence Analysis:** Odds ratios, Jaccard similarity, Fisher's exact test
- **Correlation Analysis:** Spearman ρ and Pearson *r* for pathway burden scores
- **Quadrant Classification:** Biological categorization of pathway relationships
- **Cross-Cluster Comparison:** Kruskal-Wallis tests with effect sizes (ε², Cohen's h)

### Phase 2 Key Findings

#### 1. Dose-Dependent Relationship
The **Mitochondrial-Excitotoxicity axis** emerged as the dominant mechanistic link:
- Co-occurrence: 86.5%
- Correlation: *r* = 0.687
- Odds Ratio: 38.78 (95% CI: 33.23–45.27)

This suggests cascading molecular failure where mitochondrial dysfunction drives excitotoxicity.

#### 2. Distinct Molecular Subtypes
**Proteostasis-Vesicle Trafficking** mutual exclusivity:
- Co-occurrence: 24.4%
- Correlation: *r* = −0.408
- Odds Ratio: 0.11

Indicates these pathways represent **separate disease aetiologies** rather than co-occurring mechanisms.

#### 3. Severity-Associated Gradients
| Pathway | Trend (Mild → Severe) | Cohen's h |
|:--------|:---------------------:|:---------:|
| Proteostasis | ↑ 0% → 73% | 2.05 (Large) |
| RNA Metabolism | ↑ 0% → 27% | 1.09 (Large) |
| Vesicle Trafficking | ↓ 46% → 14% | 0.73 (Medium) |
| DNA Damage | ↓ 42% → 10% | 0.77 (Medium) |

---

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/RAIDers.git
cd RAIDers

# Install dependencies
pip install -r requirements.txt
```

### Running the Full Pipeline

```bash
# Execute Phase 1 + Phase 2 pipeline
python main.py
```

### Running Individual Components

```python
# Phase 1: Federated Clustering
from fed_kmeans_components import FederatedKMeans

clusterer = FederatedKMeans(n_clusters=5, n_populations=5)
clusters = clusterer.fit(patient_features)

# Phase 2: Pathway Analysis
from phase2_pathway_annotation import PathwayAnnotator

annotator = PathwayAnnotator(gene_pathway_map='data/gene_pathway_map.json')
pathway_scores = annotator.annotate_patients(patient_variants)
```

### Interactive Dashboard

Open `als_pathway_dashboard_clusters.html` in a web browser to explore:
- Pathway prevalence by cluster
- Co-occurrence heatmaps
- Correlation matrices
- Quadrant classification scatter plot
- Network visualization

---

## 📁 Repository Structure

```
RAIDers/
├── main.py                          # Main pipeline orchestrator
├── fed_kmeans_components.py         # Federated K-means implementation
├── phase2_pathway_annotation.py     # Pathway annotation and analysis
├── stage4_pathway_analysis.py       # Statistical analysis suite
├── generate_figures.py              # Publication figure generation
├── als_pathway_dashboard.html       # Interactive visualization
├── data/
│   ├── clinvar.cleaned.csv          # Curated pathogenic variants
│   ├── gene_pathway_map.json        # Gene-to-pathway mappings
│   └── outputs/                     # Analysis results (CSV, XLSX)
├── notebooks/
│   ├── total_pipeline.ipynb         # End-to-end tutorial
│   └── biological_validation.ipynb  # Validation analyses
└── docs/
    ├── METHODOLOGY.md               # Detailed methods
    └── figures/                      # Generated figures
```

---

## 🔮 Phase 3: Multi-Institutional Collaboration (In Progress)

Phase 3 extends RAIDers to real-world multi-institutional deployment:

### Objectives
- **Real Biobank Integration:** Validation on controlled-access datasets (UK Biobank, All of Us)
- **Swarm Learning:** Transition from centralized to fully decentralized federated learning
- **Clinical Validation:** Correlation with longitudinal patient outcomes
- **Therapeutic Targeting:** Pathway-specific drug repurposing analysis

### Collaboration Opportunities

We welcome collaborations with:
- Rare disease research consortia
- Biobank data custodians
- Pharmaceutical partners for therapeutic validation
- Computational biology groups

**Contact:** [See Contributors section below]

---

## 📈 Data Sources

| Database | Purpose | Phase |
|:--------:|:--------|:-----:|
| [**ClinVar**](https://www.ncbi.nlm.nih.gov/clinvar/) | Pathogenic variant curation | 1, 2 |
| [**gnomAD**](https://gnomad.broadinstitute.org) | Population allele frequencies | 1 |
| [**OMIM**](https://www.omim.org) | Gene-disease associations | 2 |
| [**Orphanet**](https://www.orpha.net) | Rare disease ontologies | 2, 3 |
| [**Reactome**](https://reactome.org) | Pathway definitions | 2 |

---

## 📚 Citation

If you use RAIDers in your research, please cite:

```bibtex
@software{raiders2026,
  title = {RAIDers: Federated Molecular Subtyping for Rare Disease},
  author = {Shah, Aastha and Kharbanda, Arnav and others},
  year = {2026},
  url = {https://github.com/your-org/RAIDers}
}
```

---

## 👥 Contributors

| Name | Email | ORCID | Institution |
|:-----|:------|:------|:------------|
| Aastha Shah | aasthashah.work@gmail.com | [0009-0008-7811-0177](https://orcid.org/0009-0008-7811-0177) | Queen's University Belfast |
| Arnav Kharbanda | arnavkha@andrew.cmu.edu | [0009-0007-9195-9960](https://orcid.org/0009-0007-9195-9960) | Carnegie Mellon University |
| Bill Paseman | bill@rarekidneycancer.org | [0000-0002-5020-0866](https://orcid.org/0000-0002-5020-0866) | |
| Chantera Lazard | lazard.c@northeastern.edu | [0009-0006-1367-3812](https://orcid.org/0009-0006-1367-3812) | Northeastern University |
| Jialan Ma | jialanma7@gmail.com | [0009-0007-2670-9076](https://orcid.org/0009-0007-2670-9076) | Broad Institute of MIT and Harvard |
| Kushal Koirala | kkoirala@unc.edu | [0009-0009-7935-4533](https://orcid.org/0009-0009-7935-4533) | University of North Carolina |
| Kyulin Kim | lynn.kim.24@ucl.ac.uk | [0009-0007-8976-2405](https://orcid.org/0009-0007-8976-2405) | University College London |
| Nikita Rajesh | | [0009-0009-9850-5261](https://orcid.org/0009-0009-9850-5261) | Carnegie Mellon University |
| Pu (Paul) Kao | gaopuo1234@gmail.com | [0009-0003-9047-0160](https://orcid.org/0009-0003-9047-0160) | National Taiwan University |
| Shreya Nandakumar | | [0009-0006-9230-3659](https://orcid.org/0009-0006-9230-3659) | Carnegie Mellon University |
| Vibha Acharya | via16@pitt.edu | [0000-0001-6598-0052](https://orcid.org/0000-0001-6598-0052) | University of Pittsburgh |
| William Lu | wtlu@andrew.cmu.edu | [0000-0002-2768-1489](https://orcid.org/0000-0002-2768-1489) | Carnegie Mellon University |

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

<p align="center">
  <i>Developed as part of the Rare Disease & AI Initiative</i>
</p>


