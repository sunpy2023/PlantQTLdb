# <picture><source media="(prefers-color-scheme: dark)" srcset="assets/plantqtldb-header-dark.svg"><img src="assets/plantqtldb-header-light.svg" alt="PlantQTLdb" width="240"></picture>

PlantQTLdb is a plant multi-omics knowledgebase that integrates an LLM-powered research assistant. Researchers can use [PlantQTLdb](https://plantqtldb.cn/) to explore molecular variation across plant populations and evaluate genetic, molecular, and literature evidence for candidate genes and variants.

This repository provides the workflows used to analyze PlantQTLdb's multi-omics data, the [literature relationship extraction prompt](09_literature_extraction/literature_relation_extraction_prompt.md), and the [benchmark materials](10_benchmark/).

## Analysis workflows

| Directory | Analysis |
| --- | --- |
| [01_genotype](01_genotype/) | Genotype calling from whole-genome and RNA sequencing, genotype merging, and quality control |
| [02_rna_seq](02_rna_seq/) | RNA-seq preprocessing and gene expression quantification |
| [03_alternative_splicing](03_alternative_splicing/) | Splice graph construction, quantification of percent spliced in (PSI), and event classification |
| [04_alternative_polyadenylation](04_alternative_polyadenylation/) | Quantification of distal poly(A) site usage and preparation of phenotype matrices |
| [05_molqtl_mapping](05_molqtl_mapping/) | Phenotype and covariate preparation and mapping of eQTLs, sQTLs, and apaQTLs |
| [06_gwas](06_gwas/) | Genome-wide association analysis and visualization |
| [07_integrative_analysis](07_integrative_analysis/) | EWAS, TWAS, SMR, and colocalization analyses |
| [08_functional_annotation](08_functional_annotation/) | Prediction of variant consequences with SnpEff |

## Using the workflows

Choose the directory for your analysis. Before running a shell script, replace the example paths and sample identifiers with those for your dataset and install the tools it calls. Run the scripts from the repository root so that relative paths to the shared [Python and R scripts](script/) resolve correctly.

## Literature relationship extraction

The [literature extraction prompt](09_literature_extraction/literature_relation_extraction_prompt.md) defines species-aware, gene-centered relationship extraction using DeepSeek V4 Pro.

## Benchmark materials

The [benchmark directory](10_benchmark/) contains the final question set, complete responses, scores, scoring rubric and [evaluation prompt](10_benchmark/benchmark_evaluation_prompt.md). It covers 300 questions across two species groups and six task categories, evaluated under nine system configurations, including the updated MiniMax M3 results.

| File | Contents |
| --- | --- |
| [Questions](10_benchmark/questions.tsv) | 300 questions and membership of the fixed 30-question repeat subset |
| [Scoring rubric](10_benchmark/scoring_rubric.tsv) | Five evaluation dimensions and their scoring criteria (manuscript Table S5) |
| [Primary responses](10_benchmark/primary_responses.tsv) | 2,700 complete responses for generation round 1 |
| [Primary scores](10_benchmark/primary_scores.tsv) | 2,700 primary evaluations, including five dimension scores, totals and rationales (manuscript Table S6) |
| [Repeated generation responses](10_benchmark/repeated_generation_responses.tsv) | 540 additional responses from generation rounds 2 and 3 |
| [Repeated scoring](10_benchmark/repeated_scoring.tsv) | 810 evaluations from three scoring rounds of the same 270 primary responses |
| [Repeated generation scores](10_benchmark/repeated_generation_scores.tsv) | 810 evaluations of responses from three generation rounds on the fixed subset |

See the [benchmark README](10_benchmark/README.md) for the study design, data fields and distinction between repeated scoring and repeated generation. The two repeat-score files share 270 primary evaluations and should not be pooled as independent observations.
