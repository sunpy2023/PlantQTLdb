# PlantQTLdb Assistant benchmark

This directory provides the final benchmark questions, complete responses, scores, scoring rubric and evaluation prompt supporting Figure 3 and Figures S3–S5.

## Evaluation design

### Questions and task categories

The benchmark contains 300 questions: 150 for rice and 150 for *Arabidopsis thaliana*. Each of the six task categories contains 50 questions, with 25 per species group:

- **Gene function:** biological functions and mechanisms of individual genes.
- **Gene annotation:** gene identifiers and reference-genome annotations.
- **Molecular regulation:** molecular phenotypes, including gene expression, and their genetic regulation.
- **Variant analysis:** genetic variants and evidence at associated loci.
- **Candidate gene prioritization:** selection and ranking of candidate genes using genetic and functional evidence.
- **Hypothesis evaluation:** assessment of the evidence supporting a biological explanation or hypothesis.

These category names are used consistently in the data files, manuscript and figures. A fixed subset of 30 questions was used for repeated evaluations, with 15 questions per species group and five per task category.

### Systems and resource comparison

Nine system configurations were evaluated: the full PlantQTLdb Assistant, an Assistant configuration restricted to the literature-derived knowledge graph (literature-KG-only), and seven general-purpose language models: GPT-5.6 Sol, GLM-5.3, DeepSeek V4 Pro, DeepSeek V4.1 Flash, Kimi K2.5, MiniMax M3 and Qwen-3.8-Flash.

Both Assistant configurations used DeepSeek V4.1 Flash and the same agent framework. The full configuration additionally accessed population multi-omics resources and retrieval and analysis tools. This comparison assesses the combined contribution of these additional resources and tools to answer quality. Figure 3A–B and Figures S3–S5 show the full Assistant and seven general-purpose models; Figure 3C compares the two Assistant configurations. The released data retain all nine configurations.

### Response scoring

We evaluated anonymized responses using GPT-5.6 Sol on five equally weighted dimensions: factual accuracy, evidence completeness, source traceability, hallucination control and calibration, and synthesis quality. Each dimension was scored as an integer from 0 to 20; the five scores were summed to give a total of 0–100. The six task categories describe what the questions ask, whereas the five scoring dimensions assess answer quality across all categories.

The [scoring rubric](scoring_rubric.tsv) defines each dimension and its score bands (Table S5). The [evaluation prompts](benchmark_evaluation_prompt.md) specify the scoring and evidence-verification procedures. No predefined reference answer is required.

## Final data files

| File | Records | Contents |
| --- | ---: | --- |
| [questions.tsv](questions.tsv) | 300 | Question ID, species group, task category, question text and repeat-subset flag |
| [scoring_rubric.tsv](scoring_rubric.tsv) | 25 | Dimension definitions and score bands from manuscript Table S5 |
| [benchmark_evaluation_prompt.md](benchmark_evaluation_prompt.md) | — | Scoring and evidence-verification prompts |
| [primary_responses.tsv](primary_responses.tsv) | 2,700 | Complete primary responses |
| [primary_scores.tsv](primary_scores.tsv) | 2,700 | Dimension scores, totals and rationales from manuscript Table S6, linked to response IDs |
| [repeated_generation_responses.tsv](repeated_generation_responses.tsv) | 540 | Complete responses from generation rounds 2 and 3 |
| [repeated_scoring.tsv](repeated_scoring.tsv) | 810 | Three ratings of each fixed primary response in the repeat subset |
| [repeated_generation_scores.tsv](repeated_generation_scores.tsv) | 810 | One rating per response across three generation rounds in the repeat subset |

`Question ID` links each record to the question metadata in `questions.tsv`. Use `Response ID` to link a score to its complete answer in the response files; `Evaluation ID` identifies a rating. Generation and scoring rounds are recorded separately. Answers and rationales containing tabs or newlines are quoted using standard TSV/CSV conventions; physical line counts do not equal record counts.

## Repeated evaluations

- **Repeated scoring (Figure S4):** the same 270 primary responses (30 questions × nine configurations) were rated in three scoring rounds, giving 810 evaluations.
- **Repeated generation (Figure S5):** three responses were generated per question and configuration, with one rating per response, giving 810 evaluations. The 540 additional answers from rounds 2 and 3 are stored in `repeated_generation_responses.tsv`; round-1 answers are in `primary_responses.tsv`.

Round 1 in both score files reuses the primary evaluation, so the two files share 270 evaluations. These shared records should not be counted as independent repeats.

The mean scoring range is the average of the per-response maximum-minus-minimum total score across three ratings. The mean within-question generation standard deviation is calculated by taking the sample standard deviation of the three response scores for each question and averaging those 30 values within each system. The generation score variation therefore includes evaluation variability as well as differences between generated responses.
