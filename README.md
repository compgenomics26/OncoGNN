# OncoGNN

## An Interpretable Graph Neural Network Framework for Disease-Gene Prioritization Using Heterogeneous Biological Features

OncoGNN is a graph neural network (GNN)-based framework for disease-gene prioritization that integrates heterogeneous biological information within a weighted protein-protein interaction network.

The framework was developed and evaluated across three cancer types:

- Esophageal carcinoma (ESCA)
- Liver hepatocellular carcinoma (LIHC)
- Lung adenocarcinoma (LUAD)

OncoGNN integrates:

1. HPA-derived biological annotation features
2. Cancer-specific single-cell transcriptomic features
3. Network embeddings generated using node2vec
4. Weighted protein-protein interaction information from STRING

The framework evaluates multiple GNN architectures and uses probability calibration and consensus prediction to prioritize high-confidence candidate genes.

---

## Framework Overview

OncoGNN represents genes as nodes in a weighted protein-protein interaction network. Each gene is characterized by heterogeneous biological features derived from complementary data sources.

The workflow consists of:

1. Data collection and preprocessing
2. Gene identifier harmonization
3. Disease-gene label curation from DisGeNET
4. Weighted PPIN construction using STRING
5. Network representation learning using node2vec
6. Integration of biological annotation and single-cell features
7. GNN-based disease-gene prediction
8. Feature ablation analysis
9. Probability calibration
10. Consensus-based candidate gene prioritization
11. Literature validation and functional enrichment analysis

---

## Data Sources

The study uses publicly available datasets and resources.

| Resource | Purpose |
|---|---|
| Human Protein Atlas (HPA) | Gene-level biological annotation features |
| STRING | Protein-protein interaction network |
| DisGeNET | Disease-gene associations used for model labels |
| GEO | Cancer-specific single-cell RNA-sequencing datasets |
| TCGA/GDC | Genomic alteration features for additional analysis |

### Single-cell datasets

| Cancer | GEO accession |
|---|---|
| ESCA | GSE160269 |
| LIHC | GSE149614 |
| LUAD | GSE131907 |

---

## Gene Features

### HPA-derived biological annotation features

Gene-level features were obtained from the Human Protein Atlas, including biological annotations such as:

- Protein class
- Biological process
- Molecular function
- Disease involvement
- Subcellular localization

### Single-cell features

Cancer-specific single-cell RNA-sequencing datasets were processed and cell populations were consolidated into biologically meaningful cell subtypes.

Mean gene expression was calculated for each cell subtype to generate gene × cell-type feature matrices.

### Network features

Node2vec embeddings were generated from the STRING protein-protein interaction network.

The node2vec configuration used in the study was:

| Parameter | Value |
|---|---:|
| Embedding dimension | 128 |
| Return parameter (p) | 1 |
| In-out parameter (q) | 0.5 |
| Walk length | 60 |
| Walks per node | 15 |
| Context window | 10 |

---

## Disease-Gene Labeling

Disease-associated genes were obtained from DisGeNET.

Cancer-specific Gene-Disease Association (GDA) score and Evidence Index (EI) thresholds were applied to curate positive labels.

| Cancer | GDA threshold | EI threshold |
|---|---:|---:|
| ESCA | 0.1 | 0.7 |
| LIHC | 0.6 | 0.5 |
| LUAD | 0.6 | 0.5 |

Genes satisfying the corresponding thresholds were assigned positive labels, while the remaining genes were treated as unlabeled.

Disease-gene prioritization was formulated as a binary node classification problem using positive-unlabeled (PU) learning. An equal number of unlabeled genes were randomly sampled for each training trial to construct a balanced dataset.

---

## Protein-Protein Interaction Network

The protein-protein interaction network was obtained from STRING.

Only high-confidence interactions with a combined confidence score greater than 700 were retained.

STRING confidence scores were normalized to the range [0, 1] and incorporated as edge weights during graph construction.

---

## GNN Architectures

Three graph neural network architectures were evaluated:

- **GCN**
- **SGConv**
- **Weighted GraphSAGE**

The architectures were selected to compare different neighborhood aggregation strategies under a common framework.

All architectures used the same classification head, consisting of two fully connected layers with ReLU activation and dropout. Architectural differences were therefore primarily restricted to the graph convolution/message-passing layer and architecture-specific hyperparameters.

---

## Model Training

Models were implemented using PyTorch Geometric and PyTorch Lightning.

For each cancer type:

- Data were divided into training (60%), validation (20%), and test (20%) sets.
- Each architecture was trained for **25 independent trials**.
- Different trials used different random data partitions, unlabeled-gene sampling, and parameter initializations.
- Models were trained for up to 200 epochs.
- Early stopping was applied using validation loss.
- The Adam optimizer was used for model optimization.

Performance was evaluated using:

- AUROC
- Accuracy
- Precision
- Sensitivity (Recall)
- F1-score

Mean ROC curves were calculated across the 25 independent trials.

---

## Feature Ablation

The contribution of individual feature modalities was evaluated through feature ablation.

Each feature block was individually masked during inference without retraining the model.

The following feature modalities were evaluated:

- HPA-derived biological annotation features
- Single-cell features
- Node2vec network embeddings

Mean reduction in AUROC (ΔAUROC) was used to quantify the contribution of each feature block. Confidence intervals were estimated using 5,000 bootstrap resamples.

---

## Genomic Alteration Analysis

An additional analysis evaluated the contribution of cancer-specific genomic alteration features.

Genomic features were obtained from TCGA through the Genomic Data Commons (GDC), including:

- Cancer-specific mutation frequency
- Mutation frequency across GDC cohorts
- CNV gain frequency
- CNV loss frequency
- Total mutation count

SGConv was used as a common reference architecture for this analysis.

Genomic alteration features were appended to the existing feature matrix and evaluated through feature ablation.

Because genomic alteration features may overlap with information represented in existing disease-gene association resources, they were evaluated separately and were not included in the primary candidate-gene prioritization framework.

---

## Candidate Gene Prioritization

Probability calibration was evaluated using reliability diagrams for the selected architectures.

A prediction probability threshold of **0.90** was used to identify high-confidence predictions.

Candidate genes were subsequently prioritized using a calibration-guided consensus strategy across independently trained architectures.

The final candidate sets were:

| Cancer | High-confidence candidates |
|---|---:|
| ESCA | 21 |
| LIHC | 31 |
| LUAD | 24 |

Seven genes were shared across all three cancer types:

**BRCA1, HSP90AA1, HSP90AB1, NFKB1, MAPK3, SRC, and CASP3**

---

## Literature Validation

Prioritized candidate genes were systematically evaluated for supporting literature.

Gene symbols and curated HGNC aliases were combined with cancer-specific disease terms to retrieve relevant PubMed publications.

The resulting literature evidence for each candidate gene is provided in the supplementary results associated with the study.

---

## Functional Enrichment

Functional enrichment analysis was performed independently for the candidate genes identified in ESCA, LIHC, and LUAD.

The analysis included:

- KEGG pathway enrichment
- Gene Ontology (GO) Biological Process enrichment

Gene symbols were converted to Entrez Gene identifiers using `org.Hs.eg.db`.

KEGG and GO analyses were performed using `clusterProfiler`.
