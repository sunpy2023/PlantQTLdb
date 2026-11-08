# PlantQTLdb

PlantQTLdb is a plant multi-omics knowledgebase that integrates an LLM-powered research assistant. Researchers can use [PlantQTLdb](https://plantqtldb.cn/) to explore molecular variation across plant populations and evaluate genetic, molecular, and literature evidence for candidate genes and variants.

This repository currently provides the workflows used to analyze PlantQTLdb's multi-omics data and the [prompt used to extract gene-centered relationships from the literature](prompts/literature_relation_extraction_prompt.md).

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
