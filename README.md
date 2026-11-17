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

See [DEPENDENCIES.md](DEPENDENCIES.md) for the external programs, R packages, Python dependencies and input preparation. The numbered workflows are Linux/HPC templates that require the corresponding datasets and reference files.

## Literature relationship extraction

The [literature extraction prompt](09_literature_extraction/literature_relation_extraction_prompt.md) defines species-aware, gene-centered relationship extraction using DeepSeek V4 Pro. The earlier [prompts/ path](prompts/literature_relation_extraction_prompt.md) is retained for compatibility with existing links.

## Benchmark materials

The [benchmark evaluation prompt](10_benchmark/benchmark_evaluation_prompt.md) contains the complete scoring and evidence-verification instructions used in the evaluation. Literature relationship extraction uses the separate [literature extraction prompt](09_literature_extraction/literature_relation_extraction_prompt.md).

The benchmark comprises 300 questions across two species and six task categories, evaluated under nine system configurations. Manuscript Table S5 provides the scoring rubric, and Table S6 contains all 2,700 primary evaluations, including questions, five dimension scores, total scores, and scoring rationales.

| File | Contents |
| --- | --- |
| [Primary responses](10_benchmark/responses/primary_responses.tsv) | 2,700 complete responses for generation round 1 |
| [Repeated generation responses](10_benchmark/responses/repeated_generation_responses.tsv) | 540 complete responses from generation rounds 2 and 3 for the fixed 30-question subset |
| [Repeated scoring](10_benchmark/scores/repeated_scoring.tsv) | 810 evaluations: three scoring rounds of the same 270 primary responses |
| [Repeated generation scores](10_benchmark/scores/repeated_generation_scores.tsv) | 810 evaluations: one evaluation of each response from three generation rounds for the fixed 30-question subset |

The repeatability subset contains 15 questions per species and five per task category. Round 1 in both score files reuses the corresponding primary evaluations; it is not an additional scoring run. The two score files therefore share 270 round-1 evaluations and should not be pooled as independent observations. Response IDs link each evaluation to its complete answer; generation and scoring rounds are recorded separately. GPT-5.6 Sol was used as the judge. The original answer text and scoring rationales are preserved.

The [benchmark guide](10_benchmark/README.md) provides the [300-question set](10_benchmark/questions.tsv), [complete primary scores](10_benchmark/scores/primary_scores.tsv), [machine-readable rubric](10_benchmark/scoring_rubric.tsv), [statistical summaries](10_benchmark/statistics/) and reproducible validation/analysis commands. The primary score export contains all five dimension scores and rationales from manuscript Table S6. The current snapshot includes the updated MiniMax M3 answers and evaluations.

```bash
python3 10_benchmark/scripts/validate_benchmark.py
python3 -m pip install -r 10_benchmark/requirements.txt
python3 10_benchmark/scripts/analyze_benchmark.py --output-dir 10_benchmark/recomputed
```
