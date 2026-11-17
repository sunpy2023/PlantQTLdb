# Workflow dependencies and execution

The numbered directories contain Bash workflow templates for Linux/HPC environments. They are not a single automatic installation or an end-to-end workflow manager. Set the input paths, sample IDs, chromosome count, reference files and output paths at the top of the relevant script before running it. Run commands from the repository root so `script/` references resolve correctly.

## Python helpers

The analysis helpers use Python 3 and, depending on the script, NumPy, pandas, SciPy, matplotlib and pysam. Install them in an isolated environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

`requirements.txt` lists the Python dependencies; external bioinformatics programs and R packages require separate installation. The benchmark analysis has its own tested dependency pins in [10_benchmark/requirements.txt](10_benchmark/requirements.txt), while benchmark record validation needs no third-party packages.

## External programs

These versions are recorded in the workflow scripts; this table is a dependency guide rather than a lockfile for every historical dataset.

| Workflow | Programs/packages used by the templates |
| --- | --- |
| Genotype calling and QC | fastp 1.0.1, BWA-MEM2 2.2.1, HISAT2 2.2.1, SAMtools/BCFtools 1.21, Picard 3.4.0, GATK 4.6.2.0, Beagle 5.0, PLINK 1.9 and 2.0, Java |
| RNA-seq and expression | fastp, HISAT2, SAMtools, StringTie 3.0.0 |
| Alternative splicing | MAJIQ/VOILA 3.0.22 |
| Alternative polyadenylation | BEDTools 2.30.0, SAMtools, DaPars2 2.1 |
| Molecular QTLs | R, MatrixEQTL 2.3, data.table, PEER 2.0, PLINK 1.9, BCFtools |
| GWAS | PLINK 1.9, GEMMA 0.98.5 |
| EWAS | OSCA 1.20 |
| TWAS | FUSION and its R dependencies, PLINK 1.9, GCTA 1.95.1 |
| SMR | SMR 1.4.1, PLINK 1.9 |
| Colocalization | R, coloc 5.2.3, data.table, PLINK 1.9 |
| Variant annotation | SnpEff 5.4c, BCFtools, Java |

Install the software needed for the selected analysis and check that the corresponding executables are available on `PATH`; jar files and FUSION paths are configured explicitly in the templates. MAJIQ/VOILA, PEER, reference annotations and prediction resources may require their own installation or access procedures.

## Input and run preparation

Use a consistent reference assembly, chromosome naming scheme and sample identifiers across genotype, molecular phenotype, covariate and trait files. The helpers check sample alignment and input validity, but those checks do not establish biological compatibility of data from different assemblies or populations. Replace all `/path/to/` placeholders. Keep generated matrices, intermediate files and logs outside the source directories and retain the edited run configuration with each analysis.

The templates cover the analyses listed in the root README. Running the complete population pipelines requires the original sequencing, genotype, phenotype and reference inputs; those large datasets are not bundled here. The literature extraction and benchmark prompts are published separately from the model execution infrastructure.
