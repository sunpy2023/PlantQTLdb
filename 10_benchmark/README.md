# PlantQTLdb Assistant benchmark

This directory provides the final benchmark questions, complete responses, scores, scoring rubric and evaluation prompt supporting **Figure 3** and **Figures S3–S5**.

## Evaluation design

### Questions and task categories

The benchmark contains **300 questions**: **150 for rice** and **150 for *Arabidopsis thaliana***. Each of the six task categories contains **50 questions**, with **25 per species group**:

- **Gene function:** biological functions and mechanisms of individual genes.
- **Gene annotation:** gene identifiers and reference-genome annotations.
- **Molecular regulation:** molecular phenotypes, including gene expression, and their genetic regulation.
- **Variant analysis:** genetic variants and evidence at associated loci.
- **Candidate gene prioritization:** selection and ranking of candidate genes using genetic and functional evidence.
- **Hypothesis evaluation:** assessment of the evidence supporting a biological explanation or hypothesis.

These category names are used consistently in the data files, manuscript and figures. A fixed subset of **30 questions** was used for repeated evaluations, with 15 questions per species group and five per task category.

### Systems and resource comparison

**Ten system configurations** were evaluated in the primary benchmark:

- **Full PlantQTLdb Assistant**
- **Literature-KG-only Assistant:** restricted to the literature-derived knowledge graph.
- **Local KG Assistant:** accesses the local literature and multi-omics knowledge graph and PlantQTLdb's local retrieval tools. External retrieval and analysis tools are disabled for this configuration.
- **Seven general-purpose language models:** **GPT-5.6 Sol**, **GLM-5.3**, **DeepSeek V4 Pro**, **DeepSeek V4.1 Flash**, **Kimi K2.5**, **MiniMax M3** and **Qwen-3.8-Flash**.

The three Assistant configurations used **DeepSeek V4.1 Flash**. **Figure 3C** compares these configurations with DeepSeek V4.1 Flash evaluated on its own, using the same 300 questions. **Figure 3A–B** and **Figure S3** show the full Assistant and seven general-purpose models. Figures S4 and S5 include all ten configurations.

The Assistant software was updated during collection of the literature and multi-omics KG responses. Figure 3C therefore compares configurations as evaluated. The resource comparison reports means and matched-question differences with percentile 95% confidence intervals from 10,000 bootstrap samples (seed 20261116), stratified by species group and task category. Each resampled question contributes scores for all compared configurations.

### Response scoring

We evaluated anonymized responses using **GPT-5.6 Sol** on **five equally weighted dimensions**:

1. **Factual accuracy**
2. **Evidence completeness**
3. **Source traceability**
4. **Hallucination control and calibration**
5. **Synthesis quality**

Each dimension was scored as an integer from **0 to 20**. The five scores were summed to give a **total score of 0–100**. The six task categories describe what the questions ask, whereas the five scoring dimensions assess answer quality across all categories.

The [scoring rubric](scoring_rubric.tsv) defines each dimension and its score bands (Table S5). The [evaluation prompts](benchmark_evaluation_prompt.md) specify the scoring and evidence-verification procedures. No predefined reference answer is required.

## Final data files

| File | Records | Contents |
| --- | ---: | --- |
| [questions.tsv](questions.tsv) | 300 | Question ID, species group, task category, question text and repeat-subset flag |
| [scoring_rubric.tsv](scoring_rubric.tsv) | 25 | Dimension definitions and score bands from manuscript Table S5 |
| [benchmark_evaluation_prompt.md](benchmark_evaluation_prompt.md) | — | Scoring and evidence-verification prompts |
| [primary_responses.tsv](primary_responses.tsv) | 3,000 | Complete primary responses |
| [primary_scores.tsv](primary_scores.tsv) | 3,000 | Dimension scores, totals and rationales from manuscript Table S6, linked to response IDs |
| [resource_comparison.tsv](resource_comparison.tsv) | 7 | Overall and six task-category rows with means, confidence intervals and paired differences for Figure 3C |
| [repeated_generation_responses.tsv](repeated_generation_responses.tsv) | 600 | Complete responses from generation rounds 2 and 3 |
| [repeated_scoring.tsv](repeated_scoring.tsv) | 900 | Three ratings of each fixed primary response in the repeat subset |
| [repeated_generation_scores.tsv](repeated_generation_scores.tsv) | 900 | One rating per response across three generation rounds in the repeat subset |

### Linking the data files

- **`Question ID`** links each record to the question metadata in `questions.tsv`.
- **`Response ID`** links a score to its complete answer in the response files.
- **`Evaluation ID`** identifies an individual rating.

Generation and scoring rounds are recorded separately. Answers and rationales containing tabs or newlines are quoted using standard TSV/CSV conventions; physical line counts do not equal record counts.

## Repeated evaluations

Figures S4 and S5 and the three repeat files cover all three Assistant configurations and the seven general-purpose models. Both analyses use the same 30 questions.

- **Repeated scoring (Figure S4):** the same 300 primary responses (30 questions × ten configurations) were rated in three scoring rounds, giving 900 evaluations.
- **Repeated generation (Figure S5):** three responses were generated per question and configuration, with one rating per response, giving 900 evaluations. The 600 additional answers from rounds 2 and 3 are stored in `repeated_generation_responses.tsv`; round-1 answers are in `primary_responses.tsv`.

Round 1 in both score files reuses the primary evaluation, so the two files share **300 evaluations**. These shared records should not be counted as independent repeats.

### Summary measures

- **Mean scoring range:** the average of the per-response maximum-minus-minimum total score across three ratings.
- **Mean within-question generation standard deviation:** the sample standard deviation of the three response scores for each question, averaged across the 30 questions within each system.

Generation score variation includes both generation and scoring variability. For the literature and multi-omics KG configuration, the Assistant software was updated between the first and second generation rounds. Table S7 summarizes the three rounds and the within-item dispersion for each configuration.
