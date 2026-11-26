# <picture><source media="(prefers-color-scheme: dark)" srcset="assets/plantqtldb-header-dark.svg"><img src="assets/plantqtldb-header-light.svg" alt="PlantQTLdb" width="240"></picture>

[PlantQTLdb](https://plantqtldb.cn/) integrates population multi-omics data, genetic associations and literature-derived knowledge to support research on gene function and complex plant traits. **PlantQTLdb Assistant** retrieves and synthesizes evidence from these resources to answer plant research questions.

This repository provides the multi-omics analysis workflows, literature relation extraction prompt and final Assistant benchmark data accompanying PlantQTLdb.

## Repository contents

| Directory | Contents |
| --- | --- |
| [01_genotype](01_genotype/) | Genotype calling from whole-genome and RNA sequencing, genotype merging and quality control |
| [02_rna_seq](02_rna_seq/) | RNA-seq preprocessing and gene expression quantification |
| [03_alternative_splicing](03_alternative_splicing/) | Splice graph construction, percent spliced in (PSI) quantification and event classification |
| [04_alternative_polyadenylation](04_alternative_polyadenylation/) | Distal poly(A) site usage quantification and phenotype matrix preparation |
| [05_molqtl_mapping](05_molqtl_mapping/) | Phenotype and covariate preparation; eQTL, sQTL and apaQTL mapping |
| [06_gwas](06_gwas/) | Genome-wide association analysis and visualization |
| [07_integrative_analysis](07_integrative_analysis/) | EWAS, TWAS, SMR and colocalization analyses |
| [08_functional_annotation](08_functional_annotation/) | Variant consequence prediction using SnpEff |
| [09_literature_extraction](09_literature_extraction/) | Prompt for extracting species-aware, gene-centered relationships from the literature |
| [10_benchmark](10_benchmark/) | Final benchmark questions, complete responses, scores, scoring rubric and evaluation prompt |
| [script](script/) | Python and R helpers called by the analysis workflows |

## Using the analysis workflows

The shell workflows are intended for Linux/HPC environments. Choose the relevant analysis, install the tools it calls, and replace the example paths, sample identifiers and other settings with those for your dataset. Run the workflows from the **repository root** so that references to `script/` resolve correctly. Input datasets and reference files are specified in the individual workflows.

## Literature relation extraction

The [extraction prompt](09_literature_extraction/literature_relation_extraction_prompt.md) defines the relationship categories and evidence requirements used to extract gene-centered knowledge with **DeepSeek V4 Pro**.

## Assistant benchmark

- **Question set:** **300 questions**, with 150 each for rice and *Arabidopsis thaliana* and 50 in each of six task categories.
- Systems: three Assistant resource configurations and seven general-purpose language models. The Assistant configurations have access to the literature KG only, local literature and multi-omics resources, or full resources.
- **Evaluation:** responses were scored using **GPT-5.6 Sol** across five equally weighted dimensions: factual accuracy, evidence completeness, source traceability, hallucination control and calibration, and synthesis quality. Each dimension was scored from **0 to 20**, giving a **total score of 0–100**.

| File | Contents |
| --- | --- |
| [Question set](10_benchmark/questions.tsv) | 300 questions, species groups, task categories and repeat-subset membership |
| [Scoring rubric](10_benchmark/scoring_rubric.tsv) | Dimension definitions and scoring criteria from manuscript Table S5 |
| [Evaluation prompt](10_benchmark/benchmark_evaluation_prompt.md) | Instructions for scoring anonymized responses and verifying evidence |
| [Primary responses](10_benchmark/primary_responses.tsv) | 3,000 complete responses: 300 questions × ten configurations |
| [Primary scores](10_benchmark/primary_scores.tsv) | 3,000 evaluations with dimension scores, totals and rationales from manuscript Table S6 |
| [Resource comparison](10_benchmark/resource_comparison.tsv) | Overall and task-level means, 95% confidence intervals, and paired score differences for Figure 3C |
| [Repeated generation responses](10_benchmark/repeated_generation_responses.tsv) | 600 additional responses from generation rounds 2 and 3 on a fixed 30-question subset |
| [Repeated scoring](10_benchmark/repeated_scoring.tsv) | 900 evaluations: three scoring rounds of the same 300 primary responses |
| [Repeated generation scores](10_benchmark/repeated_generation_scores.tsv) | 900 evaluations: one rating per response across three generation rounds on the fixed subset |

Round 1 in both repeat-score files reuses the primary evaluation. The two files share **300 records** from the same 30 questions and ten configurations. See the [benchmark guide](10_benchmark/README.md) for model configurations, task categories, data fields and evaluation details.
