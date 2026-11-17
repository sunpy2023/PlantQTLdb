# PlantQTLdb Assistant benchmark

This directory provides the final benchmark questions, complete responses, scores, scoring rubric and evaluation prompt supporting Figure 3 and Figures S3–S5.

## Evaluation design

The benchmark contains 300 questions: 150 for rice and 150 for *Arabidopsis thaliana*. The six task categories are Gene function, Gene annotation, Molecular regulation, Variant analysis, Candidate gene prioritization, and Hypothesis evaluation. Each category contains 50 questions, including 25 per species group. The fixed repeat subset contains 30 questions, with 15 per species group and five per task category.

Nine configurations were evaluated: the full PlantQTLdb Assistant, its literature-KG-only configuration, and seven standalone models (GPT-5.6 Sol, GLM-5.3, DeepSeek V4 Pro, DeepSeek V4.1 Flash, Kimi K2.5, MiniMax M3 and Qwen-3.8-Flash). Both Assistant configurations used DeepSeek V4.1 Flash and the same agent framework. The full configuration additionally accessed population multi-omics resources and retrieval and analysis tools. Figure 3A–B and Figures S3–S5 show the full Assistant and seven standalone models; Figure 3C compares the two Assistant configurations.

GPT-5.6 Sol scored anonymized responses on five equally weighted dimensions: factual accuracy, evidence completeness, source traceability, hallucination control and calibration, and synthesis quality. Each dimension receives an integer score of 0–20; their sum is the total score of 0–100. The rubric does not require a predefined reference answer.

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

`Question ID` links questions across systems. `Response ID` identifies a generated answer; `Evaluation ID` identifies a rating. Generation and scoring rounds are recorded separately. Answers and rationales containing tabs or newlines are quoted using standard TSV/CSV conventions; physical line counts do not equal record counts. The current data include the updated MiniMax M3 primary answers and evaluations.

## Repeated evaluations

Repeated scoring holds each response fixed across three scoring rounds (Figure S4). Repeated generation evaluates three generated responses per question, with one rating per response (Figure S5). Round 1 in both files reuses the primary evaluation, so the two files share 270 records.

The mean scoring range is the average of the per-response maximum-minus-minimum total score across three ratings. The mean within-question generation standard deviation is calculated by taking the sample standard deviation of the three response scores for each question and averaging those 30 values within each system. The generation score variation therefore includes evaluation variability as well as differences between generated responses.
