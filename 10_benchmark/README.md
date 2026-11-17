# PlantQTLdb Assistant benchmark

This directory contains the question set, complete retained responses, evaluation prompt, scoring rubric, primary scores and repeat-experiment records supporting Figure 3 and Figures S3–S5.

## Study design

The primary benchmark contains 300 questions: 150 for rice and 150 for *Arabidopsis thaliana*, with 50 questions in each of six task categories. Each species/task combination contains 25 questions. The full category names in the data map to the shorter figure labels in [task_categories.tsv](task_categories.tsv).

Nine configurations were evaluated: PlantQTLdb Assistant with full resources, its literature-KG-only configuration, and seven standalone language models. Both Assistant configurations used DeepSeek V4.1 Flash and the same agent framework. The full configuration additionally accessed population multi-omics resources and retrieval and analysis tools. The seven standalone models were GPT-5.6 Sol, GLM-5.3, DeepSeek V4 Pro, DeepSeek V4.1 Flash, Kimi K2.5, MiniMax M3 and Qwen-3.8-Flash. Figure 3A–B and Figures S3–S5 show the full Assistant and seven standalone models; Figure 3C compares the two Assistant configurations.

GPT-5.6 Sol evaluated anonymized responses across five equally weighted dimensions: factual accuracy, evidence completeness, source traceability, hallucination control and calibration, and synthesis quality. Each dimension receives an integer score of 0–20; their sum is the total score of 0–100. The rubric does not require a predefined reference answer.

## Published files

| File | Records | Contents |
| --- | ---: | --- |
| [questions.tsv](questions.tsv) | 300 | One row per question, species group, category, question text and repeat-subset flag |
| [task_categories.tsv](task_categories.tsv) | 6 | Full task names and figure labels |
| [scoring_rubric.tsv](scoring_rubric.tsv) | 25 | Definitions and five score bands for each evaluation dimension; exported from manuscript Table S5 |
| [benchmark_evaluation_prompt.md](benchmark_evaluation_prompt.md) | — | Anonymous scoring and evidence-verification prompts |
| [responses/primary_responses.tsv](responses/primary_responses.tsv) | 2,700 | Complete generation-round-1 responses: 300 questions × nine configurations |
| [scores/primary_scores.tsv](scores/primary_scores.tsv) | 2,700 | Primary dimension scores, totals and rationales; exported from manuscript Table S6 and linked to response IDs |
| [responses/repeated_generation_responses.tsv](responses/repeated_generation_responses.tsv) | 540 | Additional responses for generation rounds 2 and 3: 30 questions × nine configurations × two rounds |
| [scores/repeated_scoring.tsv](scores/repeated_scoring.tsv) | 810 | Three scoring rounds for the same 270 primary responses |
| [scores/repeated_generation_scores.tsv](scores/repeated_generation_scores.tsv) | 810 | One evaluation per response from three generation rounds on the fixed 30-question subset |
| [statistics/](statistics/) | — | Recomputed means, paired differences, bootstrap confidence intervals, repeat summaries and ranks |
| [data_manifest.json](data_manifest.json) | — | Snapshot provenance, record counts and SHA-256 checksums for input files |

The primary score export preserves all five scores and rationales from Table S6. The question text and metadata are linked through `Question ID`; `System configuration` identifies the model or Assistant configuration. `Response ID` identifies one generated answer, and `Evaluation ID` identifies a rating of that answer. Generation and scoring rounds are separate fields. Tabs and newlines inside answers or rationales are quoted using standard TSV/CSV conventions; use a TSV parser rather than counting physical lines.

## Repeated scoring and generation

The fixed subset comprises 30 questions, with 15 per species group and five per task category. Its membership is recorded in `questions.tsv`.

For repeated scoring (Figure S4), each primary response is held fixed and scored in three rounds. For repeated generation (Figure S5), each system produces three responses per question, and each response is evaluated once. Round 1 in both score files reuses the primary evaluation. Consequently, the two files share 270 evaluations and must not be pooled as independent records.

The repeated-scoring summary first calculates each fixed response's maximum minus minimum total score over three ratings, then averages those ranges within each system. The repeated-generation summary first calculates the **sample standard deviation (`ddof=1`) of the three total scores for each question**, then averages the 30 standard deviations within each system. A range across baseline models, such as 6.45–9.23, is the minimum and maximum of seven model-specific mean standard deviations.

The published records are the retained evaluation outputs, including reviewed primary scores. The November 2026 snapshot incorporates the updated 300 MiniMax M3 primary answers and evaluations, its repeated-scoring evaluations, and the corresponding generation-round-1 scores. Additional generation-round-2 and -3 responses remain as retained from their original runs. The summaries analyze these recorded outputs; they do not make new model requests.

## Validate and reproduce statistics

Validation uses only the Python standard library:

```bash
python3 10_benchmark/scripts/validate_benchmark.py
python3 -m unittest discover -s 10_benchmark/tests -v
```

Recompute the statistics using Python 3.11 or later and the analysis dependencies:

```bash
python3 -m pip install -r 10_benchmark/requirements.txt
python3 10_benchmark/scripts/analyze_benchmark.py --output-dir 10_benchmark/recomputed
```

The committed statistics were generated with NumPy 2.3.5 and pandas 3.0.1; exact runtime versions and input checksums are recorded in [statistics/analysis_metadata.json](statistics/analysis_metadata.json).

The analysis uses 10,000 bootstrap resamples, sampling questions with replacement within species/task strata and preserving score pairing across systems. Confidence intervals are the 2.5th and 97.5th percentiles. Overall/baseline comparisons use seed `20261115`. Figure 3C uses seed `20261116` and task-first resampling order to reproduce the existing resource-ablation intervals. The script records both seeds and the resample count. It does not split the 300 questions into artificial replicate groups.

Output mapping:

| Manuscript item | Statistics file |
| --- | --- |
| Figure 3A | `statistics/overall_scores.tsv` (exclude literature-KG-only configuration) |
| Figure 3B | `statistics/task_scores.tsv` (exclude literature-KG-only configuration) |
| Figure 3C | `statistics/resource_ablation.tsv` |
| Figure S3 | `statistics/dimension_scores.tsv` (exclude literature-KG-only configuration) |
| Overall paired differences in Results | `statistics/paired_baseline_differences.tsv` |
| Figure S4 | `statistics/repeated_scoring_summary.tsv` and `statistics/repeated_scoring_plotted_round_means.tsv` |
| Figure S5 | `statistics/repeated_generation_scores_summary.tsv` and `statistics/repeated_generation_scores_plotted_round_means.tsv` |

Round-mean outputs include explicit ranks. Files ending in `plotted_round_means.tsv` exclude the literature-KG-only configuration and give ranks among the eight plotted systems. The corresponding `round_means.tsv` files retain all nine configurations.

## Evaluation reference

Zheng et al. (2023), [Judging LLM-as-a-Judge with MT-Bench and Chatbot Arena](https://arxiv.org/abs/2306.05685), supports the LLM-as-a-judge approach. The five equally weighted dimensions and their score bands are the PlantQTLdb benchmark rubric, rather than a scoring system prescribed by that reference.
